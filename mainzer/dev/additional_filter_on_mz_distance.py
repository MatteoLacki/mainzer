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
from math import ceil, log10

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

# settings.settings["chimeric_regression_fits_cnt"] = 0

proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df, full_matchmaker, protein_ions_matchmaker = \
        run_lipido(mz=mz,
                   intensity=intensity,
                   base_proteins=base_proteins,
                   base_lipids=base_lipids,
                   params=settings.settings,
                   verbose=verbose,
                   debug=True)

names = set(full_matchmaker.ions.query("charge == 3 and contains_protein").name)

def remove_count(name):
    split = name.split(" ")
    return split[len(split)-1]

def get_ions_without_single_bound_species(names):
    names = set(names)
    res = []
    for name in names:
        if name.count(" + ") > 1:
            protein_lipids = name.split(" + ")
            protein = protein_lipids[0]
            lipids = protein_lipids[1:]
            if any(f"{protein} + {remove_count(lipid)}" not in names for lipid in lipids):
                res.append(name)
    return res

get_ions_without_single_bound_species(names)



# plot_fit(centroids_df)

# proteins.to_csv("extended_proteins_report.csv", index=False)
# silly_prot = proteins.query("name=='Vamp2 + DOPC + 3 DOPE + RhodaminePE'")
# formula, charge = silly_prot["formula"].iloc[0], silly_prot["charge"].iloc[0]
# full_matchmaker.plot_ion_assignment(formula, charge, mz, intensity)
# full_matchmaker.ions

from mainzer.formulas import add_formulas

free_lipids = free_lipids_no_charge_df
promissing_proteins = promissing_protein_ions_df
merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
                                          promissing_protein_ions_df)