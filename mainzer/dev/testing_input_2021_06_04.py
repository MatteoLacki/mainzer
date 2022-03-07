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
import mainzer.plot
import mainzer.settings
import mainzer.lipido

data_folder = Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/small_spectrum")
# data_folder = Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/anirrudha_2021_08_05")
data_folder.exists()

spectrum_path = next(data_folder.glob("*.mzML"))
mz, intensity = read_spectrum(str(spectrum_path))
# mainzer.plot.plot_spectrum(mz, intensity)

base_lipids = read_base_lipids(data_folder/"base_lipids.csv")
base_proteins = read_base_proteins(data_folder/"base_proteins.csv")
settings = mainzer.settings.Settings.FromTOML(data_folder/"config.mainzer")
verbose = settings.settings["verbose"]

params=settings.settings
verbose=True

proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = mainzer.lipido.run_lipido(
    mz=mz,
   intensity=intensity,
   base_proteins=base_proteins,
   base_lipids=base_lipids,
   params=settings.settings,
   verbose=True)

protein_ions_matchmaker.ions.to_csv(data_folder/"proteins_only_search.csv")
ions = protein_ions_matchmaker.ions.query('single_precursor_unmatched_estimated_intensity > 10000')
plot_fit(centroids_df)


# protein_ions_matchmaker.plot_ion_assignment(
#     ion.formula.iloc[0],
#     ion.charge.iloc[0],
#     mz,
#     intensity)

# analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')
# output_folder = data_folder/spectrum_path.stem
# final_folder = output_folder/analysis_time
# final_folder.mkdir(parents=True, exist_ok=True)

# proteins.to_csv(final_folder/"proteins.csv", index=False)
# free_lipid_clusters.to_csv(final_folder/"free_lipid_clusters.csv", index=False)
# centroids_df.to_csv(final_folder/"centroids.csv")
# settings.save_toml(final_folder/"config.mainzer")

# # plot_fit(centroids_df)