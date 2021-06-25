%load_ext autoreload
%autoreload 2

from datetime import datetime
from pprint import pprint
import pandas as pd
pd.options.display.max_columns = None
import numpy as np
from pathlib import Path

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
import tqdm

data_folder = Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/small_spectrum")
data_folder.exists()
spectrum_path = next(data_folder.glob("*.mzML"))
mz, intensity = read_spectrum(str(spectrum_path))
# plot_spectrum(mz, intensity)

base_lipids = read_base_lipids(data_folder/"base_lipids.csv")
base_proteins = read_base_proteins(data_folder/"base_proteins.csv")
settings = Settings.FromTOML(data_folder/"config.mainzer")
verbose = settings.settings["verbose"]
for key,val in settings.settings.items():
    exec(key + '=val')


# proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = run_lipido(
#     mz=mz,
#    intensity=intensity,
#    base_proteins=base_proteins,
#    base_lipids=base_lipids,
#    **settings.settings)


ions = protein_ions
for key,val in regression_kwds.items():
    exec(key + '=val')
min_charge_sequence_length=min_charge_sequence_length


matchmaker = Matchmaker(ions,
                        centroids,
                        isotopic_calculator)
matchmaker.get_isotopic_summaries()
matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
matchmaker.ions
self = matchmaker

centroid_intensities, centroids_indices = self.centroids.integrated_intensity.values, self.centroids.index.values

ion_formula, ion_charge = next(iter(iter_ions_promissing))

max_expected_ppm_distance = 15

# matchmaker.assign_isotopologues_to_centroids(min_neighbourhood_intensity, verbose)
# matchmaker.estimate_max_ion_intensity(underfitting_quantile)

# after filtering if there is nothing left one has to break the calculations.