%load_ext autoreload
%autoreload 2

from datetime import datetime
from pprint import pprint
import pandas as pd
pd.options.display.max_columns = None
import numpy as np
from pathlib import Path
import tqdm
import matplotlib.pyplot as plt
import functools
import networkx as nx

from mainzer.data_frame_ops import round_df
from mainzer.isotope_ops import IsotopicCalculator
from mainzer.molecule_ops import \
    mered_proteins,\
    mered_lipids,\
    molecules2df,\
    crosslink,\
    merge_free_lipids_and_promissing_proteins
from mainzer.read import \
    read_spectrum,\
    read_base_lipids,\
    read_base_proteins
from mainzer.regression import \
    single_precursor_regression,\
    turn_single_precursor_regression_chimeric
from mainzer.signal_ops import cluster_spectrum
from mainzer.matchmaker import Matchmaker
from mainzer.settings import Settings
from mainzer.plot import plot_spectrum, plot_fit
from mainzer.lipido import run_lipido
from mainzer.matchmaker import mark_rows
from mainzer.charge_ops import cluster_charges
from mainzer.models import fit_DeconvolvedUnderfit
from mainzer.graph_ops import bipartiteGraph2regressionProblem

data_folder = Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/small_spectrum")
data_folder.exists()
spectrum_path = next(data_folder.glob("*.mzML"))
mz, intensity = read_spectrum(str(spectrum_path))
# plot_spectrum(mz, intensity)

base_lipids = read_base_lipids(data_folder/"base_lipids.csv")
base_proteins = read_base_proteins(data_folder/"base_proteins.csv")
settings = Settings.FromTOML(data_folder/"config.mainzer")
params = settings.settings
verbose = params["verbose"]


# proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = \
#         run_lipido(mz=mz,
#                    intensity=intensity,
#                    base_proteins=base_proteins,
#                    base_lipids=base_lipids,
#                    params=settings.settings,
#                    verbose=verbose)
# plot_fit(centroids_df)


#TODO: reintroduce1 Michals baseline correction idea
centroids_df = cluster_spectrum(mz, intensity)

# this is simple: filerting on `max_intensity` in `centroids_df`
filtered_centroids_df = centroids_df[
    (centroids_df.highest_intensity >= params["min_highest_intensity"]).values &\
    (centroids_df.left_mz           >= params["min_mz"]).values &\
    (centroids_df.right_mz          <= params["max_mz"]).values
].copy()

if verbose:
    print("Getting proteins mers")#TODO: change for logger
# initially we search for proteins only

protein_mers = mered_proteins(base_proteins, params["only_heteromers"])
protein_ions = molecules2df(protein_mers,
                            range(params["min_protein_cluster_charge"],
                                  params["max_protein_cluster_charge"]+1))

if verbose:
    print("Setting up IsoSpec (soon to defeat NeutronStar, if it already does not).")
isotopic_calculator = IsotopicCalculator(params["isotopic_coverage"], 
                                         params["isotopic_bin_size"])

if verbose:
    print("Checking for promissing proteins")

regression_kwds = params.copy()
regression_kwds["centroids"] = filtered_centroids_df
regression_kwds["isotopic_calculator"] = isotopic_calculator
regression_kwds["verbose"] = verbose
del regression_kwds["min_charge_sequence_length"]

protein_ions_matchmaker = \
    single_precursor_regression(protein_ions,
                                min_charge_sequence_length=params["min_charge_sequence_length"],
                                **regression_kwds)
promissing_protein_ions_df = protein_ions_matchmaker.ions[protein_ions_matchmaker.ions.reason_for_filtering_out=="none"].copy()

if verbose:
    print("Getting lipid mers")

free_lipid_mers = dict(mered_lipids(base_lipids,
                                    params["min_lipid_mers"],
                                    params["max_lipid_mers"]))
free_lipids_no_charge_df = molecules2df(free_lipid_mers)
free_lipids_with_charge_df = molecules2df(
    free_lipid_mers, 
    charges=range(params["min_free_lipid_cluster_charge"],
                  params["max_free_lipid_cluster_charge"])
)

if verbose:
    print("Getting promissing protein-lipid centroids.")


# from mainzer.formulas import add_formulas
promissing_protein_lipid_complexes = \
    merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
                                              promissing_protein_ions_df)

promissing_proteins = promissing_protein_ions_df[["name","formula","charge"]].copy()
promissing_proteins['contains_protein'] = True
promissing_protein_lipid_complexes['contains_protein'] = True
free_lipids_with_charge_df['contains_protein'] = False

# Promissing Ions = Promissing Proteins + Protein Lipid Clusters + Free Lipid Clusters
promissing_ions = pd.concat([promissing_proteins,
                             promissing_protein_lipid_complexes,
                             free_lipids_with_charge_df],
                             ignore_index=True)

full_matchmaker = \
    single_precursor_regression(promissing_ions,
                                min_charge_sequence_length=1,
                                **regression_kwds)

run_chimeric_regression = params["chimeric_regression_fits_cnt"] > 0

# full_matchmaker.ions[full_matchmaker.ions.reason_for_filtering_out == 'none']

if run_chimeric_regression:
    if verbose:
        print("Performing chimeric regression.")
    full_matchmaker = turn_single_precursor_regression_chimeric(
        full_matchmaker,
        params["fitting_to_void_penalty"], 
        merge_zeros=True,
        normalize_X=False,
        chimeric_regression_fits_cnt=params["chimeric_regression_fits_cnt"],
        min_chimeric_intensity_threshold=params["min_chimeric_intensity_threshold"],
        verbose=verbose
    )

