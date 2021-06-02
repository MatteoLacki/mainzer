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
from mainzer.lipido import run_lipido
from mainzer.plot import plot_fit, plot_spectrum
from mainzer.matchmaker import Matchmaker

from collections import defaultdict
import itertools
import json
import tqdm

# fixing clustering
with open("mainzer/dev/settings.json") as f:
    settings = json.load(f)
mz, intensity = read_spectrum(settings['path_spectrum'])


min_protein_mers = 1
max_protein_mers = 1
min_lipid_mers = 1
max_lipid_mers = 5
min_protein_cluster_charge = 1
max_protein_cluster_charge = 50
min_free_lipid_cluster_charge = 1
max_free_lipid_cluster_charge = 10
min_charge_sequence_length = 1
min_highest_intensity = 1000
min_mz = 100
max_mz = 7000
isotopic_coverage = .99
isotopic_bin_size = .1
neighbourhood_thr = 1.1
underfitting_quantile = 0.00
only_heteromers = False
min_max_intensity_threshold = 1000
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


proteins, free_lipid_clusters, centroids_df = run_lipido(# spectrum preprocessing
               mz,
               intensity,
               min_highest_intensity,
               min_mz,
               max_mz,
               # proteins
               base_proteins,
               min_protein_cluster_charge,
               max_protein_cluster_charge,
               # lipids
               base_lipids,
               min_lipid_mers,
               max_lipid_mers,
               min_free_lipid_cluster_charge,
               max_free_lipid_cluster_charge,
               # crosslinking
               only_heteromers, # discard homomers
               # filters
               min_neighbourhood_intensity,
               min_charge_sequence_length,
               min_total_fitted_probability,
               # single molecule regression
               isotopic_coverage,
               isotopic_bin_size,
               neighbourhood_thr,
               underfitting_quantile,
               min_max_intensity_threshold,
               # multiple molecule regression
               deconvolve,
               fitting_to_void_penalty,
               # verbosity might be removed in favor of a logger
               verbose=True
)

proteins, free_lipid_clusters, centroids_df
plot_fit(centroids_df)

proteins.to_csv("stringently_filtered_proteins.csv")
free_lipid_clusters.to_csv("stringently_filtered_free_lipid_clusters.csv")

test_ions = pd.read_csv("test_data/ions.csv")
matchmaker = single_precursor_regression(
    test_ions,
    min_charge_sequence_length=min_charge_sequence_length,
    **regression_kwds
)

full_matchmaker = turn_single_precursor_regression_chimeric(
    matchmaker,
    fitting_to_void_penalty, 
    merge_zeros=True,
    normalize_X=False,
    verbose=verbose
)

centroids_df2 = pd.merge(
    centroids_df,
    full_matchmaker.centroids[["chimeric_intensity_in_centroid", "chimeric_remainder"]],
    left_index=True,
    right_index=True,
    how='left'
)
plot_fit(centroids_df)

# oions = pd.read_csv("/home/matteo/Projects/och_Kallol/testres/analysis_09_04_2021__14_04_46/ions.csv")
oions = pd.read_csv("/home/matteo/Projects/och_Kallol/unlipid/testres/analysis_07_04_2021__20_30_20/ions.csv")
# olc = pd.read_csv("/home/matteo/Projects/och_Kallol/mainzer/testres/lipido__date_2021_05_13_time_15_37_23/lipid_clusters.csv")
# op = pd.read_csv("/home/matteo/Projects/och_Kallol/mainzer/testres/lipido__date_2021_05_13_time_15_37_23/proteins.csv")
# olc.sort_values("deconvolved_intensity", ascending=False).head(3)
# olc["Nycodenz_2"]
nyco2 = oions.query("name == 'Nycodenz_Nycodenz' and charge == 1")
# isotopic_calculator.spectrum("C38H52I6N6O18")

oions.query("min_isospec_mz >= 1540 and max_isospec_mz <= 1600").sort_values("deconvolved_intensity", ascending=False)

region_around_1640 = filtered_centroids_df.query("mz_apex > 1642 and mz_apex < 1650")
region_around_1640.integrated_intensity.sum()
# the intensity should amount to 2.15 million.
plot_spectrum(filtered_centroids_df.mz_apex, filtered_centroids_df.integrated_intensity)


# looks like a Nycondenz_Nycodenz^1+
free_lipids_with_charge_df.query("name == 'DOPC_2'")
len(promissing_protein_lipid_complexes)
len(free_lipids_no_charge_df) * len(promissing_protein_ions_df)

promissing_ions.query("name == 'Nycodenz_2'")

plot_spectrum(mz, intensity)
plot_spectrum(mz[(mz >= 1642.8) & (mz <= 1644.81)], intensity[(mz >= 1642.8) & (mz <= 1644.81)])

filtered_centroids_df.query("mz_apex > 1642 and mz_apex < 1645")

matchmaker = Matchmaker(promissing_ions,
                        filtered_centroids_df,
                        isotopic_calculator)
matchmaker.ions.query("name == 'Nycodenz_2' and charge == 1")

matchmaker.get_isotopic_summaries()
matchmaker.ions.query("name == 'Nycodenz_2' and charge == 1")

matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
self = matchmaker
quick=True


it = iter(promissing_ions.groupby("formula"))
formula, formula_charge = next(it)
formula_charge

%%time
for formula, formula_charge in promissing_ions.groupby("formula"):
    masses, probs = self.isotopic_calculator.masses_probs(formula)
    peaks_cnt = len(masses)
    charges_cnt = len(formula_charge)
    ion_mzs = np.tile(masses, charges_cnt)/np.repeat(formula_charge.charge.values, peaks_cnt) + self.isotopic_calculator.PROTON_MASS
    ion_probs = np.tile(probs, charges_cnt)

isotopologues2centroids.query("formula == @nyco2.formula[0] and charge == 1")
region_around_1640   
self.ion2centroids.query("formula == @nyco2.formula[0] and charge == 1")
matchmaker.assign_isotopologues_to_centroids(min_neighbourhood_intensity, verbose)
matchmaker.ions.query("name == 'Nycodenz_2' and charge == 1")
matchmaker.ion2centroids.query("formula == @nyco2.formula[0]")
matchmaker.estimate_max_ion_intensity(underfitting_quantile)
matchmaker.ions.query("name == 'Nycodenz_2' and charge == 1")
# it looks that maximal intensity estimates look right.
matchmaker.ions.query("name == 'Nycodenz_2' and charge == 1").maximal_intensity_estimate.iloc[0]
nyco2.maximal_intensity.iloc[0]



full_matchmaker.ions.query("name == 'Nycodenz_2' and charge==1")
final_ions.query("name == 'Nycodenz_2' and charge==1")
full_matchmaker

# single precursor total error is weird...
# single precursor fit error bizarre

# the maximal estimate more than 2 higher than it used to be..

# how big should it be???
centroids
plot_fit(centroids_df)







matchmaker.summarize_ion_assignments()
matchmaker.add_ion_assignments_summary_to_ions()
matchmaker.get_theoretical_intensities_where_no_signal_was_matched()
matchmaker.get_total_intensity_errors_for_maximal_intensity_estimates()
matchmaker.filter_ions_that_do_not_fit_centroids(min_total_fitted_probability)
matchmaker.filter_ions_with_low_maximal_intensity(min_max_intensity_threshold)
matchmaker.charge_ions_that_are_not_in_chargestate_sequence(min_charge_sequence_length)

full_matchmaker.ions