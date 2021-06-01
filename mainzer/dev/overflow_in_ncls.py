%load_ext autoreload
%autoreload 2

from pprint import pprint
import pandas as pd
pd.options.display.max_columns = None
import numpy as np
import aa2atom
from mainzer.molecule_ops import iter_mered_molecules, crosslink, molecules2df, iter_mers, iter_mers2, iter_mered_molecules2, mered_ions,  mered_lipids, mered_proteins, merge_free_lipids_and_promissing_proteins
from mainzer.formulas import formula2counter, counter2formula, aa2formula
from mainzer.read import read_spectrum, read_base_lipids, read_base_proteins
from mainzer.signal_ops import cluster_spectrum
from mainzer.regression import single_precursor_regression, turn_single_precursor_regression_chimeric
from mainzer.molecule_filters import charge_sequence_filter
from mainzer.isotope_ops import IsotopicCalculator
from mainzer.models import fit_DeconvolvedUnderfit
from mainzer.graph_ops import RegressionGraph


from collections import defaultdict
import itertools
import json

# fixing clustering
with open("mainzer/dev/settings.json") as f:
    settings = json.load(f)
mz, intensity = read_spectrum(settings['path_spectrum'])


min_protein_mers = 1
max_protein_mers = 3
min_lipid_mers = 1
max_lipid_mers = 4
min_protein_cluster_charge = 1
max_protein_cluster_charge = 100
min_free_lipid_cluster_charge = 1
max_free_lipid_cluster_charge = 5
min_charge_sequence_length = 3
min_highest_intensity = 100
min_mz = 100
max_mz = 18000
isotopic_coverage = .95
isotopic_bin_size = .1
neighbourhood_thr = 1.1
underfitting_quantile = 0.00
only_heteromers = False
min_max_intensity_threshold = 200
min_neighbourhood_intensity = 0
min_total_fitted_probability = .8
deconvolve = True
fitting_to_void_penalty = 1.0
verbose = True
path = "/home/matteo/Projects/och_Kallol/mainzer/test_data/molecules2.csv"
# base_lipids = read_base_lipids("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_lipids_test.csv")
# base_proteins = read_base_proteins("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_proteins_test.csv")

base_lipids = read_base_lipids("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_lipids.csv")
base_proteins = read_base_proteins("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_proteins.csv")


# clusters_df = cluster_spectrum(mz, intensity)
# filtered_clusters_df = clusters_df[(clusters_df.I_max >= min_intensity_threshold).values &\
#                                    (clusters_df.left_mz >= min_mz).values &\
#                                    (clusters_df.right_mz <= max_mz).values].copy()
# protein_mers = mered_proteins(base_proteins, only_heteromers)
# protein_ions = molecules2df(protein_mers,
#                             range(min_protein_cluster_charge,
#                                   max_protein_cluster_charge+1))
# isotopic_calculator = IsotopicCalculator(isotopic_coverage, isotopic_bin_size)

# matchmaker = single_molecule_regression(protein_ions,
#                                         filtered_clusters_df,
#                                         isotopic_calculator,
#                                         neighbourhood_thr,
#                                         underfitting_quantile,
#                                         min_total_fitted_probability,
#                                         min_max_intensity_threshold,
#                                         min_charge_sequence_length)

# promissing_protein_ions_df = matchmaker.ions
# free_lipid_mers = dict(mered_lipids(base_lipids,
#                                     min_lipid_mers,
#                                     max_lipid_mers))
# free_lipids_no_charge_df = molecules2df(free_lipid_mers)
# free_lipids_with_charge_df = \
#     molecules2df(free_lipid_mers, 
#                  charges = range(min_free_lipid_cluster_charge,
#                                  max_free_lipid_cluster_charge))
# promissing_protein_lipid_complexes = \
#     merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
#                                               promissing_protein_ions_df)
# promissing_ions = pd.concat([protein_ions,
#                              promissing_protein_lipid_complexes,
#                              free_lipids_with_charge_df],
#                              ignore_index=True)

# full_matchmaker = single_molecule_regression(promissing_ions,
#                                              filtered_clusters_df,
#                                              isotopic_calculator,
#                                              neighbourhood_thr,
#                                              underfitting_quantile,
#                                              min_total_fitted_probability,
#                                              min_max_intensity_threshold,
#                                              min_charge_sequence_length)

# # need to cut the graph to the proper things.
# full_matchmaker.build_regression_bigraph()
# full_matchmaker.fit_multiple_ion_regressions()



# full_matchmaker.promissing_ions2centroids.query("formula=='H2159C1292O360N313S10P2' and charge==6")
# full_matchmaker.promissing_ions2centroids.query("formula=='H2144C1289O361N313S10P3' and charge==6")



# this is so cool: we don't have to refit the model for non-chimeric candidates
# that saves potentially loads of time.
# full_matchmaker.promissing_ions2centroids

# nonchimeric_precursors = [chimeric_ion_group.pop() for chimeric_ion_group in full_matchmaker.G.iter_convoluted_ions() if len(chimeric_ion_group) == 1]
# nonchimeric_precursors = pd.DataFrame(nonchimeric_precursors, columns=full_matchmaker.ION)

# X = pd.merge(
#     nonchimeric_precursors,
#     full_matchmaker.promissing_ions2centroids,
#     on=full_matchmaker.ION
# )


# # check if all of things with small difference correspond to non-chimeric precursors
# full_matchmaker.promissing_ions2centroids['estim_diff'] = \
#     full_matchmaker.promissing_ions2centroids.maximal_intensity - \
#     full_matchmaker.promissing_ions2centroids.chimeric_intensity_estimate 
# full_matchmaker.promissing_ions2centroids[np.abs(full_matchmaker.promissing_ions2centroids.estim_diff) < .0000001]

# import matplotlib.pyplot as plt
# X = pd.merge(
#     full_matchmaker.chimeric_intensity_estimates,
#     full_matchmaker.max_intensity.maximal_intensity,
#     right_index=True, left_index=True
# )
# X['chimeric'] = 1
# nonchimeric_precursors_set = set(zip(nonchimeric_precursors.formula, nonchimeric_precursors.charge))
# X.chimeric.loc[[a in nonchimeric_precursors_set for a in X.index]] = 0

# plt.scatter(np.log(X.maximal_intensity),
#             np.log(X.chimeric_intensity_estimate),
#             c=X.chimeric)
# plt.show()





# matchmaker2.estimate_max_ion_intensity(underfitting_quantile)
# matchmaker2.summarize_ion_assignments()
# matchmaker2.add_ion_assignments_summary_to_ions()
# matchmaker2.get_theoretical_intensities_where_no_signal_was_matched()
# matchmaker2.get_total_intensity_errors_for_maximal_intensity_estimates()
# matchmaker2.filter_ions_that_do_not_fit_centroids(min_total_fitted_probability)
# matchmaker2.filter_ions_with_low_maximal_intensity(min_max_intensity_threshold)
# matchmaker2.charge_ions_that_are_not_in_chargestate_sequence(min_charge_sequence_length)

# this in not quick! pandas likes bigger data passed once.
# # OK, so we can do it more efficient!
# # Great: what about plotting the outcomes?

# # simple_centroids = self.centroids.I_sum.copy().reset_index()

# # intensities = self.centroids.I_sum.values
# # clusters = self.centroids.index.values
# # formula, charge = next( promissing_ions.itertuples(index=False) )
# def iter_local_ion2centroids(intensities, clusters):
#     for formula, charge in promissing_ions.itertuples(index=False):
#         mz, probs = self.isotopic_calculator.spectrum(formula, charge)
#         edges = self.centroids_intervals.point_query(mz)
#         if len(edges.query):
#             local_isotopologue2centoids = pd.DataFrame({
#                 "formula": formula,
#                 "charge": charge,
#                 "cluster": clusters[edges.interval_db],
#                 "I_sum": intensities[edges.interval_db],
#                 "mz": mz[edges.query],
#                 "prob": probs[edges.query]
#             })
#             local_isotopologue2centoids['mzprob'] = local_isotopologue2centoids.mz * local_isotopologue2centoids.prob
#             local_ion2centroids = local_isotopologue2centoids.groupby(self.ION+["cluster","I_sum"]).agg(
#                 {"mz":[min, max],
#                  "mzprob":[sum, len],
#                  "prob":sum}
#             )
#             local_ion2centroids.columns = ["min_mz", "max_mz", "tot_mzprob","isotopologues_cnt", "isospec_isotopologues_prob"]
#             local_ion2centroids['weighted_mz'] = local_ion2centroids.tot_mzprob / local_ion2centroids.isospec_isotopologues_prob
#             local_ion2centroids = local_ion2centroids.reset_index()
#             local_ion2centroids.isotopologues_cnt = local_ion2centroids.isotopologues_cnt.astype(np.uint64)
#             local_ion2centroids["intensity_over_prob"] = local_ion2centroids.I_sum / local_ion2centroids.isospec_isotopologues_prob
#             yield local_ion2centroids

# ion2centroids = pd.concat(iter_local_ion2centroids(self.centroids.I_sum.values, self.centroids.index.values))