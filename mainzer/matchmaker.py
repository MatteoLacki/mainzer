import functools
import numpy as np
import pandas as pd
import tqdm

from dataclasses import dataclass
from typing import Tuple, List

from .intervals import IntervalQuery
from .graph_ops import RegressionGraph
from .models import fit_DeconvolvedUnderfit
from .charge_ops import cluster_charges

import mainzer.isotope_ops


def mark_rows(
    table: pd.DataFrame,
    ions_to_mark,
    names,
    column_name: str,
    marking,
    alternative: bool=None
) -> pd.DataFrame:
    table = table.set_index(names)
    table.loc[ions_to_mark, column_name] = marking
    if alternative is not None:
        table[column_name].fillna(alternative, inplace=True)
    return table.reset_index()


@dataclass
class Matchmaker:
    ions: pd.DataFrame
    centroids: pd.DataFrame
    isotopic_calculator: mainzer.isotope_ops.IsotopicCalculator
    ION: Tuple[str] = ("formula","charge")

    def __post_init__(self):
        assert all(colname in self.ions.columns for colname in ("name","formula","charge")), "'ions_df' should be a pandas.DataFrame with columns 'name', 'formula', and 'charge'."
        # self.ions["reason_for_filtering_out"] = np.nan
        self.centroids_intervals = IntervalQuery(self.centroids.left_mz,
                                                 self.centroids.right_mz)
        self._ion_lst = list(self.ION)

    def get_isotopic_summaries(self) -> None:
        """These include summaries of one ion's envelope: its size, total probability, minimal and maximal m/z, top-probable m/z."""
        self.ions = self.isotopic_calculator.ions_summary(self.ions)


    def get_neighbourhood_intensities(
        self,
        neighbourhood_thr: float=1.1
    ) -> None:
        """Get neighbourhood intensities.

        Get the signal intensity above the minimal m/z of the envelope minus neighbourhood_thr and below the maximal m/z of the envelope plus neighbourhood_thr.
    
        Arguments:
            neighbourhood_thr (float): The threshold for neighbourhood.
        """
        edges = self.centroids_intervals.interval_query(
            self.ions.envelope_min_mz - neighbourhood_thr,
            self.ions.envelope_max_mz + neighbourhood_thr
        )
        # different neighbourhoods might include intensity from the same centroids.
        self.ions["neighbourhood_intensity"] = self.centroids.integrated_intensity.iloc[edges.interval_db].groupby(edges.query).sum()
        self.ions.neighbourhood_intensity.fillna(0, inplace=True)# NaN -> 0


    def assign_isotopologues_to_centroids(
        self, 
        min_neighbourhood_intensity: float=0,
        verbose: bool=False
    ) -> None:
        """Assign isotopologues to centroids."""

        # using pandas it is faster to provide larger chunks of data to C++.
        # In C++ one could efficiently implement this code totally differently: asking for presence of isotopologues one after another.
        ions_promissing = self.ions[self.ions.neighbourhood_intensity > 0][self._ion_lst]
        ions_promissing.drop_duplicates(inplace=True)

        def iter_local_isotopologues2centroids():
            iter_ions_promissing = ions_promissing.itertuples(index=False)
            if verbose:
                print("Assigning isotopologues to centroids.")
                iter_ions_promissing = tqdm.tqdm(iter_ions_promissing,
                                                 total=len(ions_promissing))
            for ion_formula, ion_charge in iter_ions_promissing:
                mz, probs = self.isotopic_calculator.spectrum(ion_formula, ion_charge)
                edges = self.centroids_intervals.point_query(mz)
                if len(edges.query):
                    yield pd.DataFrame({
                        "formula":  ion_formula,
                        "charge":   ion_charge,
                        "centroid": self.centroids.index.values[edges.interval_db],
                        "mz":       mz[edges.query],
                        "prob":     probs[edges.query]
                    })

        # alternative to make graph here immediately
        self.isotopologues2centroids = pd.concat(
            iter_local_isotopologues2centroids(),
            ignore_index=True
        )

        self.isotopologues2centroids['mz_apex'] = self.centroids.mz_apex.loc[self.isotopologues2centroids.centroid].values

        self.isotopologues2centroids['signed_absolute_distance_apex_isotopologue'] = \
            self.isotopologues2centroids.mz - self.isotopologues2centroids.mz_apex

        self.isotopologues2centroids['ppm_distance_apex_isotopologue'] = \
            1e6 * \
            self.isotopologues2centroids.signed_absolute_distance_apex_isotopologue / \
            self.isotopologues2centroids.mz_apex.abs()

        self.isotopologues2centroids['mzprob'] = \
            self.isotopologues2centroids.mz * self.isotopologues2centroids.prob
        
        self.ions2centroids = self.isotopologues2centroids.groupby(self._ion_lst + ["centroid", "mz_apex"]).agg(
            {   "mz":     [min, max],
                "mzprob": [sum, len],
                "prob":    sum
            })
        self.ions2centroids.columns = [
            "min_mz", 
            "max_mz",
            "tot_mzprob",
            "isotopologues_cnt",
            "isotopologue_prob"
        ]
        self.ions2centroids.isotopologues_cnt = self.ions2centroids.isotopologues_cnt.astype(np.int64)

        self.ions2centroids.drop(columns="tot_mzprob", inplace=True)

        self.ions2centroids = self.ions2centroids.reset_index(["centroid","mz_apex"])

        self.ions2centroids['integrated_intensity'] = self.centroids.integrated_intensity.loc[self.ions2centroids.centroid].values

        self.ions2centroids["intensity_over_prob"] = self.ions2centroids.integrated_intensity / self.ions2centroids.isotopologue_prob

        indices_of_highest_theoretical_peaks = self.isotopologues2centroids.groupby(self._ion_lst+["centroid"]).prob.idxmax().values

        self.ions2centroids["top_mz"] = \
            self.isotopologues2centroids.mz.iloc[indices_of_highest_theoretical_peaks].values

        self.ions2centroids["top_prob"] = \
            self.isotopologues2centroids.prob.iloc[indices_of_highest_theoretical_peaks].values

        self.ions2centroids.top_prob /= self.ions2centroids.groupby(self._ion_lst).top_prob.sum()

        self.ions2centroids["ppm_distance_top_prob_mz_cluster_apex"] = \
            self.isotopologues2centroids.ppm_distance_apex_isotopologue.iloc[indices_of_highest_theoretical_peaks].values

        self.ions2centroids["weighted_ppm_distance_top_prob_mz_cluster_apex"] = \
            self.ions2centroids.ppm_distance_top_prob_mz_cluster_apex *\
            self.ions2centroids.top_prob

        self.ions = mark_rows(
            self.ions,
            ions_to_mark=[ion for ion in zip(self.ions.formula, self.ions.charge) if ion not in self.ions2centroids.index],
            names=self._ion_lst,
            column_name="reason_for_filtering_out",
            marking="not_within_any_centroid",
            alternative="none"
        )
        # I dunno which key I want to have, so I just don't have any and select one when necessary, and go back to no keys when I am finished.
        self.ions2centroids = self.ions2centroids.reset_index()


    def plot_ion_assignment(
        self,
        formula: str,
        charge: int,
        mz: np.array,
        intensity: np.array,
        mz_buffer: float=50,
        show: bool=True
    ) -> None:
        """This function is used for debugging purposes."""
        ION = self.isotopologues2centroids[(self.isotopologues2centroids.formula == formula) & (self.isotopologues2centroids.charge == charge) ]
        if len(ION):            
            from .plot import plot_spectrum
            import matplotlib.patches as mpatches
            import matplotlib.pyplot as plt
            MZ = mz[ (ION.mz.min() - mz_buffer <= mz) & (ION.mz.max() + mz_buffer >= mz) ] 
            I = intensity[ (ION.mz.min() - mz_buffer <= mz) & (ION.mz.max() + mz_buffer >= mz) ]
            plot_spectrum(MZ, I, show=False)
            CENTROIDS_AROUND = self.centroids[self.centroids.mz_apex.between(ION.mz.min() - mz_buffer, ION.mz.max() + mz_buffer)]
            for centroid in CENTROIDS_AROUND.itertuples():
                rect = mpatches.Rectangle((centroid.left_mz, 0),
                                          centroid.right_mz - centroid.left_mz,
                                          centroid.highest_intensity, 
                                          fill=False,
                                          color="purple",
                                          linewidth=2)
                plt.gca().add_patch(rect)
            HIT_CENTROIDS = self.centroids.loc[ ION.centroid.unique() ]
            for centroid in HIT_CENTROIDS.itertuples():
                rect = mpatches.Rectangle((centroid.left_mz, 0),
                                          centroid.right_mz - centroid.left_mz,
                                          centroid.highest_intensity, 
                                          fill=False,
                                          color="black",
                                          linewidth=3)
                plt.gca().add_patch(rect)
            ion_mzs, ion_probs = self.isotopic_calculator.spectrum(formula, charge)
            plt.scatter(ion_mzs, ion_probs/ion_probs.max() * HIT_CENTROIDS.highest_intensity.max(), c='orange' )
            if show:
                plt.show()


    def mark_ions_far_from_apexes_of_centroids(
        self, 
        max_expected_ppm_distance: float=15
    ) -> None:
        """This function only marks thing."""
        expected_ppm_dist = self.ions2centroids.groupby(self._ion_lst).weighted_ppm_distance_top_prob_mz_cluster_apex.sum()
        ions_far_from_centroid_apexes = expected_ppm_dist[expected_ppm_dist.abs() > max_expected_ppm_distance]

        mark_far_away_ions = functools.partial(
            mark_rows, 
            ions_to_mark=ions_far_from_centroid_apexes.index,
            names=self._ion_lst,
            column_name="reason_for_filtering_out",
            marking="far_from_centroid_apex")

        self.ions2centroids = mark_far_away_ions(self.ions2centroids)
        self.isotopologues2centroids = mark_far_away_ions(self.isotopologues2centroids)
        self.ions = mark_far_away_ions(self.ions)
        self.ions.set_index(self._ion_lst, inplace=True)
        self.ions["expected_ppm_dist"] = expected_ppm_dist
        self.ions.expected_ppm_dist.fillna(-1, inplace=True)
        self.ions.reset_index(inplace=True)

    def estimate_max_ion_intensity(
        self,
        underfitting_quantile: float=0
    ) -> None:
        assert underfitting_quantile >= 0 and underfitting_quantile <= 1, "'underfitting_quantile' is a probability and so must be between 0 and 1."
        
        # self.ions = self.ions.set_index(self._ion_lst)
        # self.ions['maximal_intensity_estimate'] = self.ions2centroids.groupby(self._ion_lst).intensity_over_prob.quantile(underfitting_quantile)
        # self.ions['maximal_intensity_estimate'].fillna(0, inplace=True)


        self.maximal_intensity_estimate = self.ions2centroids.groupby(self._ion_lst).intensity_over_prob.quantile(underfitting_quantile)
        self.maximal_intensity_estimate.name = "maximal_intensity_estimate"
        self.ions2centroids = self.ions2centroids.merge(
            self.maximal_intensity_estimate,
            left_on=self._ion_lst, 
            right_index=True,
            how="left"
        )
        # getting errors:
        self.ions2centroids["single_precursor_fit_error"] = \
            self.ions2centroids.maximal_intensity_estimate * \
            self.ions2centroids.isotopologue_prob

        single_precursor_fit_error = self.ions2centroids.groupby(self._ion_lst).single_precursor_fit_error.sum()

        self.maximal_intensity_estimate = pd.concat([
            self.maximal_intensity_estimate,
            single_precursor_fit_error
        ], axis=1)
        # feed maximal_intensity and its error to ions
        self.ions = self.ions.merge(
            self.maximal_intensity_estimate,
            left_on=self._ion_lst, 
            right_index=True,
            how="left"
        )
        self.ions.fillna(0, inplace=True)


    def summarize_ion_assignments(self) -> None:
        self.ion_assignments_summary = self.ions2centroids.groupby(self._ion_lst).agg(
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
            self.ions[self._ion_lst + ["envelope_total_prob"]],
            left_index=True,
            right_on=self._ion_lst
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

        self.ions = pd.merge(
            self.ions,
            self.ion_assignments_summary.drop("envelope_total_prob", axis=1),
            on=self._ion_lst,
            how="left"
        )
        self.ions.fillna(0, inplace=True)
        self.ions.explainable_centroids = self.ions.explainable_centroids.astype(int)


    def get_theoretical_intensities_where_no_signal_was_matched(self) -> None:
        # another way to measure error: the intensity of unmatched fragments
        self.ions["single_precursor_unmatched_estimated_intensity"] =\
            self.ions.maximal_intensity_estimate * \
            self.ions.envelope_unmatched_prob


    def get_total_intensity_errors_for_maximal_intensity_estimates(self) -> None:
        self.ions["single_precursor_total_error"] = \
            self.ions.single_precursor_fit_error + \
            self.ions.single_precursor_unmatched_estimated_intensity


    def mark_ions_with_probability_mass_outside_centroids(
        self,
        min_total_fitted_probability: float=.8,
    ) -> None:
        self.ions.reason_for_filtering_out = np.where(
            (self.ions.reason_for_filtering_out == "none") & \
            (self.ions.envelope_matched_to_signal < min_total_fitted_probability), 
            "probability_mass_mostly_outside_centroids",
            self.ions.reason_for_filtering_out
        )


    def mark_ions_with_low_maximal_intensity(
        self,
        min_max_intensity_threshold: float=0,
    ) -> None:
        self.ions.reason_for_filtering_out = np.where(
            (self.ions.reason_for_filtering_out == "none") & 
            (self.ions.maximal_intensity_estimate < min_max_intensity_threshold), 
            "low_maximal_intensity_estimate",
            self.ions.reason_for_filtering_out
        )


    def mark_ions_not_in_charge_cluster(
        self,
        min_charge_sequence_length: int=1,
    ) -> None:
        if min_charge_sequence_length > 1:
            self.ions = self.ions.set_index(self._ion_lst)
            formulas_charges = self.ions[self.ions.reason_for_filtering_out == "none"].index.to_frame(False)
            formulas_charges["charge_group"] = formulas_charges.groupby("formula").charge.transform(cluster_charges)
            formulas_charges = formulas_charges.set_index(self._ion_lst)
            self.ions["charge_group"] = formulas_charges
            self.ions.charge_group.fillna(0, inplace=True)
            self.ions.charge_group = self.ions.charge_group.astype(int)
            self.ions.reason_for_filtering_out = np.where(
                (self.ions.reason_for_filtering_out == "none") & 
                (self.ions.charge_group == 0), 
                "not_in_charge_group",
                self.ions.reason_for_filtering_out
            )
            self.ions.reset_index(inplace=True)


    def build_regression_bigraph(self) -> None:
        good_ions = self.ions[self.ions.reason_for_filtering_out == 'none'][["name"]+self._ion_lst]
        self.promissing_ions2centroids = pd.merge(
            self.ions2centroids,
            good_ions,
            on=self._ion_lst
        )
        
        promissing_ions_assignments_summary = pd.merge(
            self.ion_assignments_summary,
            good_ions,
            on=self._ion_lst
        )

        self.G = RegressionGraph()
        for r in self.promissing_ions2centroids.itertuples():
            ion = r.formula, r.charge
            self.G.add_edge(ion, r.centroid, prob=r.isotopologue_prob )
            self.G.nodes[r.centroid]["intensity"] = r.integrated_intensity

        for r in promissing_ions_assignments_summary.itertuples():
            ion = r.formula, r.charge
            self.G.nodes[ion]["prob_out"] = r.envelope_unmatched_prob


    def get_chimeric_ions(self) -> List:
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
        self.chimeric_intensity_estimates.index.names = self._ion_lst

        self.ions = pd.merge(# passing estimates unto ions
            self.ions,
            self.chimeric_intensity_estimates,
            left_on=self._ion_lst,
            right_index=True,
            how="left"
        )
        self.promissing_ions2centroids = pd.merge(# passing estimates unto ions
            self.promissing_ions2centroids,
            self.chimeric_intensity_estimates,
            left_on=self._ion_lst,
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


    def filter_out_ions_with_low_chimeric_estimates(self, min_chimeric_intensity_threshold):
        self.ions = self.ions[self.ions.chimeric_intensity_estimate >= min_chimeric_intensity_threshold]


    def reset_state_to_before_chimeric_regression(self):
        """This needs to be run to repeat regression [with lower number of control variables]."""
        self.ions.drop(columns=["chimeric_intensity_estimate"],
                       inplace=True)
        self.centroids.drop(columns=["chimeric_remainder",
                                     "chimeric_intensity_in_centroid"], 
                            inplace=True)
        del self.models
        del self.chimeric_intensity_estimates
        del self.promissing_ions2centroids
        del self.G


    def get_chimeric_groups(self):
        chimeric_groups = pd.concat(self.G.iter_chimeric_groups(), ignore_index=True)
        assert len(chimeric_groups) == len(self.ions)
        self.ions = pd.merge(
            self.ions,
            chimeric_groups,
            on=self._ion_lst
        )

