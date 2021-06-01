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
from mainzer.centroiding import QueryCentroids
from mainzer.regression import single_molecule_regression, Matchmaker
from mainzer.molecule_filters import charge_sequence_filter
from mainzer.isotope_ops import IsotopicCalculator
from mainzer.graph_ops import get_regression_bigraph
from mainzer.models import fit_DeconvolvedUnderfit



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
min_intensity_threshold = 100
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


clusters_df = cluster_spectrum(mz, intensity)
filtered_clusters_df = clusters_df[(clusters_df.I_max >= min_intensity_threshold).values &\
                                   (clusters_df.left_mz >= min_mz).values &\
                                   (clusters_df.right_mz <= max_mz).values].copy()
protein_mers = mered_proteins(base_proteins, only_heteromers)
protein_ions = molecules2df(protein_mers,
                            range(min_protein_cluster_charge,
                                  max_protein_cluster_charge+1))
isotopic_calculator = IsotopicCalculator(isotopic_coverage, isotopic_bin_size)

matchmaker = single_molecule_regression(protein_ions,
                                        filtered_clusters_df,
                                        isotopic_calculator,
                                        neighbourhood_thr,
                                        underfitting_quantile,
                                        min_total_fitted_probability,
                                        min_max_intensity_threshold,
                                        min_charge_sequence_length)

promissing_protein_ions_df = matchmaker.ions
free_lipid_mers = dict(mered_lipids(base_lipids,
                                    min_lipid_mers,
                                    max_lipid_mers))
free_lipids_no_charge_df = molecules2df(free_lipid_mers)
free_lipids_with_charge_df = \
    molecules2df(free_lipid_mers, 
                 charges = range(min_free_lipid_cluster_charge,
                                 max_free_lipid_cluster_charge))
promissing_protein_lipid_complexes = \
    merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
                                              promissing_protein_ions_df)
promissing_ions = pd.concat([protein_ions,
                             promissing_protein_lipid_complexes,
                             free_lipids_with_charge_df],
                             ignore_index=True)

matchmaker2 = Matchmaker(promissing_ions, filtered_clusters_df, isotopic_calculator)
matchmaker2.get_isotopic_summaries()
matchmaker2.get_neighbourhood_intensities(neighbourhood_thr)

self = matchmaker2

promissing_ions = self.ions[self.ions.neighbourhood_intensity > 0][self.ION]
promissing_ions.drop_duplicates(inplace=True)

# MAIN IDEA: make separate queries for each isotopic envelope
# full_envelopes = pd.concat((self.isotopic_calculator.to_frame(form, z)
#                             for form, z in promissing_ions.itertuples(index=False)),
#                            ignore_index=True)

# it takes some time to query centroids


def iter_local_isotopologues2centroids(intensities, clusters):
    for formula, charge in promissing_ions.itertuples(index=False):
        mz, probs = self.isotopic_calculator.spectrum(formula, charge)
        edges = self.centroids_intervals.point_query(mz)
        if len(edges.query):
            yield pd.DataFrame({
                "formula": formula,
                "charge": charge,
                "cluster": clusters[edges.interval_db],
                "I_sum": intensities[edges.interval_db],
                "mz": mz[edges.query],
                "prob": probs[edges.query]
            })


isotopologues2centroids = pd.concat(iter_local_isotopologues2centroids(self.centroids.I_sum.values, self.centroids.index.values), ignore_index=True)


edges = self.centroids_intervals.point_query(full_envelopes.mz.values)
print(edges.query.min()) # wow
print(edges.query.max()) # wow


isotopologue2centroid = pd.concat([
    full_envelopes.iloc[edges.query].reset_index(drop=True),
    self.centroids[important_columns].iloc[edges.interval_db].reset_index()
], axis=1)





# get an work-around for this interval class...




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