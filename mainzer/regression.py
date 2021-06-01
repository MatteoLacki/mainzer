import numpy as np
import pandas as pd
import tqdm

from .intervals import IntervalQuery
from .graph_ops import RegressionGraph
from .models import fit_DeconvolvedUnderfit
from .molecule_filters import charge_sequence_filter

DEBUG = True


class Matchmaker(object):
    def __init__(self,
                 ions,
                 centroids,
                 isotopic_calculator):
        assert all(colname in ions.columns for colname in ("name","formula","charge")), "'ions_df' should be a pandas.DataFrame with columns 'name', 'formula', and 'charge'."
        self.ions = ions
        self.centroids = centroids
        self.centroids_intervals = IntervalQuery(self.centroids.left_mz,
                                                 self.centroids.right_mz)
        self.isotopic_calculator = isotopic_calculator
        self.ION = ["formula","charge"]

    def get_isotopic_summaries(self):
        self.ions = self.isotopic_calculator.ions_summary(self.ions)

    def get_neighbourhood_intensities(self, neighbourhood_thr=1.1):
        edges = self.centroids_intervals.interval_query(self.ions.isospec_min_mz + neighbourhood_thr,
                                                        self.ions.isospec_max_mz + neighbourhood_thr)
        # different neighbourhoods might include intensity from the same centroids.
        self.ions["neighbourhood_intensity"] = self.centroids.I_sum.iloc[edges.interval_db].groupby(edges.query).sum()
        self.ions.neighbourhood_intensity.fillna(0, inplace=True)# NaN = 0 

    def assign_isotopologues_to_centroids(self, min_neighbourhood_intensity=0):
        promissing_ions = self.ions[self.ions.neighbourhood_intensity > 0][self.ION]
        promissing_ions.drop_duplicates(inplace=True)

        # using pandas it is faster to provide larger chunks of data to C++.
        # In C++ one could efficiently implement this code totally differently: asking for presence of isotopologues one after another.
        def iter_local_isotopologues2centroids(intensities, centroids):
            for formula, charge in promissing_ions.itertuples(index=False):
                mz, probs = self.isotopic_calculator.spectrum(formula, charge)
                edges = self.centroids_intervals.point_query(mz)
                if len(edges.query):
                    yield pd.DataFrame({
                        "formula": formula,
                        "charge": charge,
                        "centroid": centroids[edges.interval_db],
                        "I_sum": intensities[edges.interval_db],
                        "mz": mz[edges.query],
                        "prob": probs[edges.query]
                    })
        #TODO: alternative to make graph here immediately
        isotopologues2centroids = \
            pd.concat(iter_local_isotopologues2centroids(self.centroids.I_sum.values, 
                                                         self.centroids.index.values),
                      ignore_index=True)
        isotopologues2centroids['mzprob'] = \
            isotopologues2centroids.mz * isotopologues2centroids.prob
        
        # here we add I_sum to the key because it is unique to each centroid.
        self.ion2centroids = isotopologues2centroids.groupby(self.ION+["centroid","I_sum"]).agg(
            {"mz":[min, max],
             "mzprob":[sum, len],
             "prob":sum}
        )
        self.ion2centroids.columns = ["min_mz", "max_mz", "tot_mzprob","isotopologues_cnt", "isospec_isotopologues_prob"]
        self.ion2centroids['weighted_mz'] = self.ion2centroids.tot_mzprob / self.ion2centroids.isospec_isotopologues_prob
        self.ion2centroids = self.ion2centroids.reset_index()
        self.ion2centroids.isotopologues_cnt = self.ion2centroids.isotopologues_cnt.astype(np.int64)
        self.ion2centroids["intensity_over_prob"] = self.ion2centroids.I_sum / self.ion2centroids.isospec_isotopologues_prob
        
    def estimate_max_ion_intensity(self, underfitting_quantile=0):
        assert underfitting_quantile >= 0 and underfitting_quantile <= 1, "'underfitting_quantile' is a probability and so must be between 0 and 1."
        self.max_intensity = self.ion2centroids.groupby(self.ION).intensity_over_prob.quantile(underfitting_quantile)
        self.max_intensity.name = "maximal_intensity"
        self.ion2centroids = self.ion2centroids.merge(
            self.max_intensity,
            left_on=self.ION, 
            right_index=True,
            how="left"
        )
        # getting errors:
        self.ion2centroids["max_intensity_matched_fit_error"] = \
            self.ion2centroids.maximal_intensity * \
            self.ion2centroids.isospec_isotopologues_prob

        max_intensity_matched_fit_error = self.ion2centroids.groupby(self.ION).max_intensity_matched_fit_error.sum()

        self.max_intensity = pd.concat([
            self.max_intensity,
            max_intensity_matched_fit_error
        ], axis=1)
        # feed maximal_intensity and its error to ions
        self.ions = self.ions.merge(
            self.max_intensity,
            left_on=self.ION, 
            right_index=True,
            how="left"
        )
        self.ions.fillna(0, inplace=True)

    def summarize_ion_assignments(self):
        self.ion_assignments_summary = self.ion2centroids.groupby(self.ION).agg(
            {"isospec_isotopologues_prob":"sum",
             "I_sum":"sum", # Here sum is OK: we sum different centroids' intensities.
             "centroid":"nunique"}
        )
        self.ion_assignments_summary = self.ion_assignments_summary.rename(
            columns={"isospec_isotopologues_prob":"isospec_prob_with_signal",
                     "I_sum":"proximity_intensity", 
                     "centroid":"touched_centroids"}
        )
        self.ion_assignments_summary = pd.merge(
            self.ion_assignments_summary,
            self.ions[self.ION + ["isospec_final_coverage"]],
            left_index=True,
            right_on=self.ION
        )
        self.ion_assignments_summary["isospec_prob_without_signal"] = \
            self.ion_assignments_summary.isospec_final_coverage - \
            self.ion_assignments_summary.isospec_prob_with_signal

        # sometimes the numerics kick in here and produce small negative numbers.
        isospec_prob_without_signal = self.ion_assignments_summary.isospec_prob_without_signal
        # checking if they ain't big.
        assert np.all(np.abs(isospec_prob_without_signal[isospec_prob_without_signal < 0]) < 1e-10), "The unfitted probabilties are negative and big."
        # fix by setting to 0
        self.ion_assignments_summary.isospec_prob_without_signal = np.where(isospec_prob_without_signal < 0, 0, isospec_prob_without_signal)

    def add_ion_assignments_summary_to_ions(self):
        self.ions = pd.merge(
            self.ions,
            self.ion_assignments_summary.drop("isospec_final_coverage", axis=1),
            on=self.ION,
            how="left"
        )
        self.ions.fillna(0, inplace=True)
        self.ions.touched_centroids = self.ions.touched_centroids.astype(int)

    def get_theoretical_intensities_where_no_signal_was_matched(self):
        # another way to measure error: the intensity of unmatched fragments
        self.ions["max_intensity_unmatched_fit_error"] =\
            self.ions.maximal_intensity * \
            self.ions.isospec_prob_without_signal

    def get_total_intensity_errors_for_maximal_intensity_estimates(self):
        self.ions["max_intensity_total_error"] = \
            self.ions.max_intensity_matched_fit_error + \
            self.ions.max_intensity_unmatched_fit_error

    def filter_ions_with_low_maximal_intensity(self, min_max_intensity_threshold):
        self.ions = self.ions[self.ions.maximal_intensity >= min_max_intensity_threshold]

    def filter_ions_that_do_not_fit_centroids(self, min_total_fitted_probability=.8):
        self.ions = self.ions[self.ions.isospec_prob_with_signal >= min_total_fitted_probability]

    def charge_ions_that_are_not_in_chargestate_sequence(self, min_charge_sequence_length=1):
        if min_charge_sequence_length > 1:
            self.ions = charge_sequence_filter(self.ions,
                                               min_charge_sequence_length)

    def build_regression_bigraph(self):
        self.promissing_ions2centroids = \
            pd.merge(self.ion2centroids,
                     self.ions[["name"]+self.ION],
                     on=self.ION)

        self.promissing_ions_assignments_summary = \
            pd.merge(self.ion_assignments_summary,
                     self.ions[["name"]+self.ION],
                     on=self.ION)

        self.G = RegressionGraph()
        for r in self.promissing_ions2centroids.itertuples():
            ion = r.formula, r.charge
            self.G.add_edge(ion, r.centroid, prob=r.isospec_isotopologues_prob )
            self.G.nodes[r.centroid]["intensity"] = r.I_sum

        for r in self.promissing_ions_assignments_summary.itertuples():
            ion = r.formula, r.charge
            self.G.nodes[ion]["prob_out"] = r.isospec_prob_without_signal

    def get_chimeric_ions(self):
        return list(self.G.iter_convoluted_ions())

    def fit_multiple_ion_regressions(self,
                                     fitting_to_void_penalty=1.0, 
                                     merge_zeros=True,
                                     normalize_X=False,
                                     verbose=True):
        #TODO: might not run estimates for single candidates: they are already calculated and 
        # but the problem is that the user might set different quantile.
        # so don't let him.
        iter_problems = self.G.iter_regression_problems(merge_zeros, normalize_X)
        if verbose:
            iter_problems = tqdm.tqdm(iter_problems,
                                      total=self.G.count_connected_components())
        
        self.models = [fit_DeconvolvedUnderfit(X,Y, lam=fitting_to_void_penalty)
                       for X,Y in iter_problems]

        self.chimeric_intensity_estimates = pd.concat([m.coef for m in self.models])
        self.chimeric_intensity_estimates.name = 'chimeric_intensity_estimate'
        self.chimeric_intensity_estimates.index.names = self.ION

        self.ions = pd.merge(# passing estimates unto ions
            self.ions,
            self.chimeric_intensity_estimates,
            left_on=self.ION,
            right_index=True,
            how="left"
        )
        self.promissing_ions2centroids = pd.merge(# passing estimates unto ions
            self.promissing_ions2centroids,
            self.chimeric_intensity_estimates,
            left_on=self.ION,
            right_index=True,
            how="left"
        )
        self.promissing_ions2centroids["chimeric_peak_intensity"] = \
            self.promissing_ions2centroids.chimeric_intensity_estimate * \
            self.promissing_ions2centroids.isospec_isotopologues_prob

        self.chimeric_intensity_in_centroid = self.promissing_ions2centroids.groupby("centroid").chimeric_peak_intensity.sum()
        self.chimeric_intensity_in_centroid.name = "chimeric_intensity_in_centroid"
        self.centroids = pd.merge(
            self.centroids,
            self.chimeric_intensity_in_centroid,
            left_index=True,
            right_index=True,
            how='left'
        )
        self.centroids.chimeric_intensity_in_centroid.fillna(0.0, inplace=True)
        self.centroids["chimeric_remainder"] = \
            self.centroids.I_sum - self.centroids.chimeric_intensity_in_centroid

        assert np.all(np.abs(self.centroids.chimeric_remainder[self.centroids.chimeric_remainder < 0]) < 1e-10), "The unfitted probabilties are negative and big."
        self.centroids.chimeric_remainder = np.where(
                self.centroids.chimeric_remainder < 0,
                0, 
                self.centroids.chimeric_remainder)

    @staticmethod
    def from_existing_centroids(existing_matchmaker, new_ions):
        """This is introduced to save time on unnecessary initialization."""
        return IonsCentroids(new_ions,
                             existing_matchmaker.centroids,
                             existing_matchmaker.isotopic_calculator)



def single_precursor_regression(ions,
                                centroids,
                                isotopic_calculator,
                                neighbourhood_thr=1.1,
                                underfitting_quantile=0.0,
                                min_total_fitted_probability=.8,
                                min_max_intensity_threshold=100,
                                min_charge_sequence_length=1):
    matchmaker = Matchmaker(ions, centroids, isotopic_calculator)
    matchmaker.get_isotopic_summaries()
    matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
    matchmaker.assign_isotopologues_to_centroids()
    matchmaker.estimate_max_ion_intensity(underfitting_quantile)
    matchmaker.summarize_ion_assignments()
    matchmaker.add_ion_assignments_summary_to_ions()
    matchmaker.get_theoretical_intensities_where_no_signal_was_matched()
    matchmaker.get_total_intensity_errors_for_maximal_intensity_estimates()
    matchmaker.filter_ions_that_do_not_fit_centroids(min_total_fitted_probability)
    matchmaker.filter_ions_with_low_maximal_intensity(min_max_intensity_threshold)
    matchmaker.charge_ions_that_are_not_in_chargestate_sequence(min_charge_sequence_length)
    return matchmaker


def chimeric_regression(ions,
                        centroids,
                        isotopic_calculator,
                        neighbourhood_thr=1.1,
                        underfitting_quantile=0.0,
                        min_total_fitted_probability=.8,
                        min_max_intensity_threshold=100,
                        min_charge_sequence_length=1,
                        fitting_to_void_penalty=1.0, 
                        merge_zeros=True,
                        normalize_X=False,
                        verbose=True):
    matchmaker = Matchmaker(ions, centroids, isotopic_calculator)
    matchmaker.get_isotopic_summaries()
    matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
    matchmaker.assign_isotopologues_to_centroids()
    matchmaker.estimate_max_ion_intensity(underfitting_quantile)
    matchmaker.summarize_ion_assignments()
    matchmaker.add_ion_assignments_summary_to_ions()
    matchmaker.get_theoretical_intensities_where_no_signal_was_matched()
    matchmaker.get_total_intensity_errors_for_maximal_intensity_estimates()
    matchmaker.filter_ions_that_do_not_fit_centroids(min_total_fitted_probability)
    matchmaker.filter_ions_with_low_maximal_intensity(min_max_intensity_threshold)
    matchmaker.charge_ions_that_are_not_in_chargestate_sequence(min_charge_sequence_length)
    matchmaker.build_regression_bigraph()
    matchmaker.get_chimeric_ions()
    matchmaker.fit_multiple_ion_regressions()
    return matchmaker


def turn_single_precursor_regression_chimeric(matchmaker,
                                              fitting_to_void_penalty=1.0, 
                                              merge_zeros=True,
                                              normalize_X=False,
                                              verbose=True):
    matchmaker.build_regression_bigraph()
    matchmaker.get_chimeric_ions()
    matchmaker.fit_multiple_ion_regressions()
    return matchmaker
