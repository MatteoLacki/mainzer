import numpy as np
import pandas as pd
import tqdm

from .intervals import IntervalQuery
from .graph_ops import get_regression_bigraph
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
        #TODO: one might optimize this for graph purposes
        promissing_ions = self.ions[self.ions.neighbourhood_intensity > 0][self.ION]
        promissing_ions.drop_duplicates(inplace=True)
        full_envelopes = pd.concat((self.isotopic_calculator.to_frame(form, z)
                                    for form, z in promissing_ions.itertuples(index=False)),
                                   ignore_index=True)
        
        # edges = self.centroids_intervals.point_query(full_envelopes.head(50).mz)
        # full_envelopes.iloc[edges.query]

        # i = 0 
        # for i in range(len(full_envelopes)):
        #     edges = self.centroids_intervals.point_query(full_envelopes.mz.iloc[(i*100):((i+1)*100)])
        #     full_envelopes.iloc[edges.query].reset_index(drop=True)
        # np.ulonglong
        # np.arange()

        edges = self.centroids_intervals.point_query(full_envelopes.mz.values)
        print(edges.query.min()) # wow
        print(edges.query.max()) # wow

        important_columns = ["I_sum","left_mz","right_mz"]

        isotopologue2centroid = pd.concat([
            full_envelopes.iloc[edges.query].reset_index(drop=True),
            self.centroids[important_columns].iloc[edges.interval_db].reset_index()
        ], axis=1)

        # assert all((isotopologue2centroid.mz >= isotopologue2centroid.left_mz).values & (isotopologue2centroid.mz <= isotopologue2centroid.right_mz).values ), "Some theoretical m/z values are not inside  intervals matched to centroids."
        isotopologue2centroid['mzprob'] = isotopologue2centroid.mz * isotopologue2centroid.prob
        # here we add I_sum to the key because it is unique to each cluster.
        self.ion2centroids = isotopologue2centroid.groupby(self.ION+["cluster","I_sum"]).agg(
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
             "cluster":"nunique"}
        )
        self.ion_assignments_summary = self.ion_assignments_summary.rename(
            columns={"isospec_isotopologues_prob":"isospec_prob_with_signal",
                     "I_sum":"proximity_intensity", 
                     "cluster":"touched_centroids"}
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

    @staticmethod
    def from_existing_centroids(existing_matchmaker, new_ions):
        """This is introduced to save time on unnecessary initialization."""
        return IonsCentroids(new_ions,
                             existing_matchmaker.centroids,
                             existing_matchmaker.isotopic_calculator)



def single_molecule_regression(ions,
                               clusters,
                               isotopic_calculator,
                               neighbourhood_thr=1.1,
                               underfitting_quantile=0.0,
                               min_total_fitted_probability=.8,
                               min_max_intensity_threshold=100,
                               min_charge_sequence_length=1):
    matchmaker = Matchmaker(ions, clusters, isotopic_calculator)
    matchmaker.get_isotopic_summaries()
    matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
    self = matchmaker

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



def single_molecule_regression_old(centroids,
                                   ions_df,
                                   isotopic_coverage=.95,
                                   isotopic_bin_size=.1,
                                   neighbourhood_thr=1.1,
                                   underfitting_quantile=0.05,
                                   verbose=False,
                                   **kwds):

    assert all(colname in ions_df.columns for colname in ("name","formula","charge")), "'ions_df' should be a pandas.DataFrame with columns 'name', 'formula', and 'charge'."
    assert isotopic_coverage >= 0 and isotopic_coverage <= 1, "'isotopic_coverage' is a probability and so must be between 0 and 1."
    assert underfitting_quantile >= 0 and underfitting_quantile <= 1, "'underfitting_quantile' is a probability and so must be between 0 and 1."

    if verbose:
        print("Getting Isotopic Envelopes")

    isotopic_envelopes = IsotopicEnvelopes(ions_df.formula.unique(),
                                           isotopic_coverage,
                                           isotopic_bin_size)
    if DEBUG:
        print("Ctored")

    ions_df = ions_df.merge(isotopic_envelopes.charged_envelopes_summary(ions_df.formula, ions_df.charge))
    
    if DEBUG:
        print("merged")
    
    ion_idx = ['formula','charge']
    ions_df = ions_df.set_index(ion_idx)

    if verbose:
        print(f"Getting intensity of signal centroids in a {neighbourhood_thr}-neighbourhood of theoretical peaks.")

    minmax_signals = centroids.interval_query(ions_df.min_isospec_mz - neighbourhood_thr,
                                              ions_df.max_isospec_mz + neighbourhood_thr)
    ions_df["neighbourhood_intensity"] = minmax_signals.I_sum.groupby(ion_idx).sum()
    ions_df.neighbourhood_intensity.fillna(0, inplace=True)# NaN = 0 intensity
    del minmax_signals


    if verbose:
        print("Assigning isotopic envelope peaks to real signal centroids.")

    full_envelopes = isotopic_envelopes.to_frame(ions_df[ions_df.neighbourhood_intensity > 0].index)
    peak_assignments = centroids.point_query(full_envelopes.isospec_mz)
    peak_assignments = pd.concat([full_envelopes.iloc[peak_assignments.index],
                                  peak_assignments], axis=1)
    del full_envelopes

    assert np.all(peak_assignments.isospec_mz == peak_assignments.query_mz), "Something is wrong while creating peak assignments: query m/z not the same as isospec_mz."
    peak_assignments.drop("query_mz", axis=1, inplace=True) # "query_mz" is the same as "isospec_mz"
    peak_assignments = peak_assignments[peak_assignments.I_sum > 0]
    # need to have a measure of error

    # signal-based centroiding of isospec peaks! so we do not really need that binning!!!
    # watch out for summing "I_sum" twice! "I_sum" must be part of a group here:
    peak_assignments_clustered = peak_assignments.groupby(ion_idx + ["cluster","left_mz","right_mz","mz_apex","I_sum"]).agg({"isospec_prob":"sum", "isospec_mz":['min','max']}).reset_index()
    peak_assignments_clustered.columns = [f"{a}_{b}" if a == "isospec_mz" else a for a, b in peak_assignments_clustered.columns] # how I hate pandas..


    if verbose:
        print("Getting maximal intensity estimates.")

    max_intensity = peak_assignments_clustered.set_index(ion_idx)
    max_intensity = max_intensity.I_sum / max_intensity.isospec_prob
    max_intensity = max_intensity.groupby(ion_idx).quantile(underfitting_quantile)
    ions_df["maximal_intensity"] = max_intensity
    ions_df.maximal_intensity.fillna(0, inplace=True)


    if verbose:
        print("Summarizing assignments.")

    peak_assignments_summary = peak_assignments_clustered.groupby(ion_idx).agg({"isospec_prob":"sum", "I_sum":"sum", "cluster":"nunique"}).rename(columns={"isospec_prob":"isospec_prob_with_signal", "I_sum":"proximity_intensity", "cluster":"touched_centroids"})
    peak_assignments_summary["isospec_final_coverage"] = ions_df.isospec_final_coverage
    peak_assignments_summary["isospec_prob_without_signal"] = \
        peak_assignments_summary.isospec_final_coverage - \
        peak_assignments_summary.isospec_prob_with_signal

    # adding peak assignment summary to ions_df
    ions_df = pd.merge(ions_df,
                    peak_assignments_summary[['isospec_prob_with_signal',
                                              'isospec_prob_without_signal',
                                              'proximity_intensity',
                                              'touched_centroids']], 
                    left_index=True, 
                    right_index=True,
                    how='left')
    ions_df.isospec_prob_with_signal.fillna(0, inplace=True)
    ions_df.proximity_intensity.fillna(0, inplace=True)
    ions_df.touched_centroids.fillna(0, inplace=True)
    ions_df.touched_centroids = ions_df.touched_centroids.astype(int)
    ions_df.isospec_prob_without_signal = np.where(ions_df.isospec_prob_without_signal.isna(), ions_df.isospec_final_coverage, ions_df.isospec_prob_without_signal)

    return ( ions_df, 
             centroids,
             peak_assignments_clustered,
             peak_assignments_summary )



def multiple_molecule_regression(ions_df,
                                 centroids,
                                 peak_assignments_clustered, 
                                 peak_assignments_summary,
                                 fitting_to_void_penalty=1.0,
                                 verbose=False):
    ion_idx = ['formula','charge']

    if verbose:
        print("Building deconvolution graph.")

    G = get_regression_bigraph(peak_assignments_clustered, peak_assignments_summary)
    convoluted_ions_df = list(G.iter_convoluted_ions_df()) # groups of ions_df competing for signal explanation


    if verbose:
        print("Fitting multiple molecule regressions_df.")
    iter_problems = G.iter_regression_problems(merge_zeros=True,
                                               normalize_X=False) # think about this normalization

    if verbose:
        iter_problems = tqdm.tqdm(iter_problems, total=len(convoluted_ions_df))
    
    models = [fit_DeconvolvedUnderfit(X,Y, lam=fitting_to_void_penalty)
              for X,Y in iter_problems]

    estimates = pd.concat([m.coef for m in models])
    estimates.name = 'estimate'
    estimates.index.names = ion_idx

    ions_df["deconvolved_intensity"] = estimates
    ions_df.deconvolved_intensity.fillna(0, inplace=True)

    # getting errors
    peak_assignments_clustered = peak_assignments_clustered.merge(estimates,
                                                                  left_on=ion_idx,
                                                                  right_index=True).rename(columns={"estimate":"alpha"})
    peak_assignments_clustered.alpha.fillna(0, inplace=True)
    peak_assignments_clustered["under_estimate"] = \
        peak_assignments_clustered.alpha * peak_assignments_clustered.isospec_prob
    
    centroids.df["under_estimate"] = peak_assignments_clustered.groupby("cluster").under_estimate.sum()
    centroids.df.under_estimate.fillna(0.0, inplace=True)
    centroids.df["under_estimate_remainder"] = centroids.df.I_sum - centroids.df.under_estimate

    # assert np.all(ions_df.maximal_intensity >= ions_df.deconvolved_intensity + eps), "Maximal estimates lower than deconvolved!"
    #TODO: check if the difference is small: eps

    if verbose:
        print("Multiple regression done!")
    return ions_df, centroids, peak_assignments_clustered, G, convoluted_ions_df, models
