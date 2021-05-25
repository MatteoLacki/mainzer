# Use this to debug run_lipido
%load_ext autoreload
%autoreload 2

import json
import pandas as pd
import numpy as np
import tqdm

from mainzer.read import read

with open("mainzer/dev/settings.json") as f:
    settings = json.load(f)

mz, intensity = read(settings['path_spectrum'])
molecules = pd.read_csv(settings['path_molecules'])

min_protein_mers = 2
max_protein_mers = 2
max_lipid_mers = 4
min_protein_charge = 1
max_protein_charge = 10
min_lipid_charge = 1
max_lipid_charge = 10

isotopic_coverage = .95
isotopic_bin_size = .1
neighbourhood_thr = 1.1
underfitting_quantile = 0.05

deconvolve = True
fitting_to_void_penalty = 1.0

verbose = True

from mainzer.ion_generators import get_lipido_ions
from mainzer.centroiding import centroid_spectrum
from mainzer.deconv import single_molecule_regression, multiple_molecule_regression
pd.options.display.max_columns = None

ions.to_csv("test_ions.csv")

# debug get_lipid_ions
import aa2atom
from mainzer.molecule_ops import iter_mered_molecules, crosslink, molecules2df

molecules = pd.read_csv("/home/matteo/Projects/och_Kallol/mainzer/test_data/molecules_big.csv")
allowed_column_names = "group", "name", "sequence_or_formula", "min_mers", "max_mers"
assert all(col in molecules.columns for col in allowed_column_names), f"The csv with molecules should have columns: {', '.join(allowed_column_names)}."
proteins = molecules[molecules.group == "protein"]
protein_formulas = dict(zip(proteins.name, 
                            [aa2atom.aa2atom(fasta) for fasta in proteins.sequence_or_formula]))


lipids = molecules.query("group == 'lipid'")
lipids = dict(zip(lipids.name, lipids.sequence_or_formula))
lipid_mers = dict(iter_mered_molecules(lipids, max_lipid_mers))
lipid_protein_mers = dict(crosslink(lipid_mers, proteins))



