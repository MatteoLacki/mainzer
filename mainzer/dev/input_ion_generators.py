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
from mainzer.regression import single_molecule_regression
from mainzer.molecule_filters import charge_sequence_filter
from mainzer.isotope_ops import IsotopicCalculator
from mainzer.graph_ops import get_regression_bigraph
from mainzer.models import fit_DeconvolvedUnderfit
from mainzer.intervals import IntervalQuery, IonsCentroids


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
minimal_maximal_intensity_threshold = 200
min_neighbourhood_intensity = 0

deconvolve = True
fitting_to_void_penalty = 1.0

verbose = True
path = "/home/matteo/Projects/och_Kallol/mainzer/test_data/molecules2.csv"
# base_lipids = read_base_lipids("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_lipids_test.csv")
# base_proteins = read_base_proteins("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_proteins_test.csv")

base_lipids = read_base_lipids("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_lipids.csv")
base_proteins = read_base_proteins("/home/matteo/Projects/och_Kallol/mainzer/test_data/base_proteins.csv")


