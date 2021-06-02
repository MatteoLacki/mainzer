import numpy as np
import pandas as pd
import tqdm

from .intervals import IntervalQuery
from .graph_ops import RegressionGraph
from .models import fit_DeconvolvedUnderfit
from .molecule_filters import charge_sequence_filter

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

    def get_neighbourhood_intensities(self,
                                      neighbourhood_thr=1.1):
        edges = self.centroids_intervals.interval_query(
            self.ions.envelope_min_mz - neighbourhood_thr,
            self.ions.envelope_max_mz + neighbourhood_thr
        )
        # different neighbourhoods might include intensity from the same centroids.
        self.ions["neighbourhood_intensity"] = self.centroids.integrated_intensity.iloc[edges.interval_db].groupby(edges.query).sum()
        self.ions.neighbourhood_intensity.fillna(0, inplace=True)# NaN -> 0

    def assign_isotopologues_to_centroids(self, 
                                          min_neighbourhood_intensity=0,
                                          verbose=False,
                                          quick=True):
        ions_promissing = self.ions[self.ions.neighbourhood_intensity > 0][self.ION]
        ions_promissing.drop_duplicates(inplace=True)

        # if quick:
        #     self.isotopic_calculator
        #     pass
        # else:
        # using pandas it is faster to provide larger chunks of data to C++.
        # In C++ one could efficiently implement this code totally differently: asking for presence of isotopologues one after another.
        def iter_local_isotopologues2centroids(intensities, centroids):
            iter_ions_promissing = ions_promissing.itertuples(index=False)
            if verbose:
                print("Assigning isotopologues to centroids.")
                iter_ions_promissing = tqdm.tqdm(iter_ions_promissing,
                                                 total=len(ions_promissing))
            for formula, charge in iter_ions_promissing:
                mz, probs = self.isotopic_calculator.spectrum(formula, charge)
                edges = self.centroids_intervals.point_query(mz)
                if len(edges.query):
                    yield pd.DataFrame({
                        "formula": formula,
                        "charge": charge,
                        "centroid": centroids[edges.interval_db],
                        "integrated_intensity": intensities[edges.interval_db],
                        "mz": mz[edges.query],
                        "prob": probs[edges.query]
                    })
        #TODO: alternative to make graph here immediately
        isotopologues2centroids = \
            pd.concat(iter_local_isotopologues2centroids(self.centroids.integrated_intensity.values, 
                                                         self.centroids.index.values),
                      ignore_index=True)

        isotopologues2centroids['mzprob'] = \
            isotopologues2centroids.mz * isotopologues2centroids.prob
        
        # here we add integrated_intensity to the key because it is unique to each centroid.
        self.ion2centroids = isotopologues2centroids.groupby(self.ION+["centroid","integrated_intensity"]).agg(
            {"mz":[min, max],
             "mzprob":[sum, len],
             "prob":sum}
        )
        self.ion2centroids.columns = ["min_mz", "max_mz", "tot_mzprob","isotopologues_cnt", "isotopologue_prob"]
        self.ion2centroids['weighted_mz'] = self.ion2centroids.tot_mzprob / self.ion2centroids.isotopologue_prob
        self.ion2centroids.drop(columns="tot_mzprob", inplace=True)
        self.ion2centroids = self.ion2centroids.reset_index()
        self.ion2centroids.isotopologues_cnt = self.ion2centroids.isotopologues_cnt.astype(np.int64)
        self.ion2centroids["intensity_over_prob"] = self.ion2centroids.integrated_intensity / self.ion2centroids.isotopologue_prob
        
    def estimate_max_ion_intensity(self,
                                   underfitting_quantile=0):
        assert underfitting_quantile >= 0 and underfitting_quantile <= 1, "'underfitting_quantile' is a probability and so must be between 0 and 1."
        self.maximal_intensity_estimate = self.ion2centroids.groupby(self.ION).intensity_over_prob.quantile(underfitting_quantile)
        self.maximal_intensity_estimate.name = "maximal_intensity_estimate"
        self.ion2centroids = self.ion2centroids.merge(
            self.maximal_intensity_estimate,
            left_on=self.ION, 
            right_index=True,
            how="left"
        )
        # getting errors:
        self.ion2centroids["single_precursor_fit_error"] = \
            self.ion2centroids.maximal_intensity_estimate * \
            self.ion2centroids.isotopologue_prob

        single_precursor_fit_error = self.ion2centroids.groupby(self.ION).single_precursor_fit_error.sum()

        self.maximal_intensity_estimate = pd.concat([
            self.maximal_intensity_estimate,
            single_precursor_fit_error
        ], axis=1)
        # feed maximal_intensity and its error to ions
        self.ions = self.ions.merge(
            self.maximal_intensity_estimate,
            left_on=self.ION, 
            right_index=True,
            how="left"
        )
        self.ions.fillna(0, inplace=True)

    def summarize_ion_assignments(self):
        self.ion_assignments_summary = self.ion2centroids.groupby(self.ION).agg(
            {"isotopologue_prob":"sum",
             "integrated_intensity":"sum", # Here sum is OK: we sum different centroids' intensities.
             "centroid":"nunique"}
        )
        self.ion_assignments_summary = self.ion_assignments_summary.rename(
            columns={"isotopologue_prob":"envelope_matched_to_signal",
                     "integrated_intensity":"envelope_proximity_intensity", 
                     "centroid":"explainable_centroids"}
        )
        self.ion_assignments_summary = pd.merge(
            self.ion_assignments_summary,
            self.ions[self.ION + ["envelope_total_prob"]],
            left_index=True,
            right_on=self.ION
        )
        self.ion_assignments_summary["envelope_unmatched_prob"] = \
            self.ion_assignments_summary.envelope_total_prob - \
            self.ion_assignments_summary.envelope_matched_to_signal

        # sometimes the numerics kick in here and produce small negative numbers.
        envelope_unmatched_prob = self.ion_assignments_summary.envelope_unmatched_prob
        # checking if they ain't big.
        assert np.all(np.abs(envelope_unmatched_prob[envelope_unmatched_prob < 0]) < 1e-10), "The unfitted probabilties are negative and big."
        # fix by setting to 0
        self.ion_assignments_summary.envelope_unmatched_prob = np.where(envelope_unmatched_prob < 0, 0, envelope_unmatched_prob)

    def add_ion_assignments_summary_to_ions(self):
        self.ions = pd.merge(
            self.ions,
            self.ion_assignments_summary.drop("envelope_total_prob", axis=1),
            on=self.ION,
            how="left"
        )
        self.ions.fillna(0, inplace=True)
        self.ions.explainable_centroids = self.ions.explainable_centroids.astype(int)

    def get_theoretical_intensities_where_no_signal_was_matched(self):
        # another way to measure error: the intensity of unmatched fragments
        self.ions["single_precursor_unmatched_estimated_intensity"] =\
            self.ions.maximal_intensity_estimate * \
            self.ions.envelope_unmatched_prob

    def get_total_intensity_errors_for_maximal_intensity_estimates(self):
        self.ions["single_precursor_total_error"] = \
            self.ions.single_precursor_fit_error + \
            self.ions.single_precursor_unmatched_estimated_intensity

    def filter_ions_with_low_maximal_intensity(self, min_max_intensity_threshold):
        self.ions = self.ions[self.ions.maximal_intensity_estimate >= min_max_intensity_threshold]

    def filter_ions_that_do_not_fit_centroids(self, min_total_fitted_probability=.8):
        self.ions = self.ions[self.ions.envelope_matched_to_signal >= min_total_fitted_probability]

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
            self.G.add_edge(ion, r.centroid, prob=r.isotopologue_prob )
            self.G.nodes[r.centroid]["intensity"] = r.integrated_intensity

        for r in self.promissing_ions_assignments_summary.itertuples():
            ion = r.formula, r.charge
            self.G.nodes[ion]["prob_out"] = r.envelope_unmatched_prob

    def get_chimeric_ions(self):
        return list(self.G.iter_convoluted_ions())

    def fit_multiple_ion_regressions(self,
                                     fitting_to_void_penalty=1.0, 
                                     merge_zeros=True,
                                     normalize_X=False,
                                     verbose=False):
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
            self.promissing_ions2centroids.isotopologue_prob

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
            self.centroids.integrated_intensity - self.centroids.chimeric_intensity_in_centroid

        assert np.all(np.abs(self.centroids.chimeric_remainder[self.centroids.chimeric_remainder < 0]) < 1e-4), "The unfitted probabilties are negative and big."
        self.centroids.chimeric_remainder = np.where(
                self.centroids.chimeric_remainder < 0,
                0, 
                self.centroids.chimeric_remainder)

    @staticmethod
    def from_existing_centroids(existing_matchmaker, new_ions):
        """This is introduced to save time on unnecessary initialization."""
        #TODO: stop recalculating centroid intervals.
        return IonsCentroids(new_ions,
                             existing_matchmaker.centroids,
                             existing_matchmaker.isotopic_calculator)

