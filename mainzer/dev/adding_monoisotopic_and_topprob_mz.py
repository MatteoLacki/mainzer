%load_ext autoreload
%autoreload 2

from datetime import datetime
from pprint import pprint
import pandas as pd
pd.options.display.max_columns = None
import numpy as np
from pathlib import Path


from mainzer.molecule_ops import iter_mered_molecules, crosslink, molecules2df, iter_mers, iter_mers2, iter_mered_molecules2, mered_ions,  mered_lipids, mered_proteins, merge_free_lipids_and_promissing_proteins
from mainzer.formulas import formula2counter, counter2formula, aa2formula
from mainzer.read import read_spectrum, read_base_lipids, read_base_proteins
from mainzer.settings import Settings
from mainzer.lipido import run_lipido
from mainzer.plot import plot_spectrum, plot_fit

from mainzer.signal_ops import cluster_spectrum
from mainzer.isotope_ops import IsotopicCalculator
from mainzer.molecule_ops import \
    mered_proteins, \
    mered_lipids, \
    molecules2df, \
    crosslink, \
    merge_free_lipids_and_promissing_proteins
from mainzer.read import read_spectrum, read_base_lipids, read_base_proteins
from mainzer.regression import single_precursor_regression, turn_single_precursor_regression_chimeric
from mainzer.matchmaker import Matchmaker

data_folder = Path("test_data/2021_06_04")
spectrum_path = next(data_folder.glob("*.mzML"))
mz, intensity = read_spectrum(str(spectrum_path))

settings = Settings.FromTOML(data_folder/"config.mainzer")
for key,val in settings.settings.items():
    exec(key + '=val')
# plot_spectrum(mz, intensity)

base_lipids = read_base_lipids(data_folder/"base_lipids.csv")
base_proteins = read_base_proteins(data_folder/"base_proteins.csv")




output_folder = Path("test_data/2021_06_04")