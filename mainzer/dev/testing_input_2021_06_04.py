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

data_folder = Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/small_spectrum")
data_folder.exists()
spectrum_path = next(data_folder.glob("*.mzML"))
mz, intensity = read_spectrum(str(spectrum_path))
# plot_spectrum(mz, intensity)

base_lipids = read_base_lipids(data_folder/"base_lipids.csv")
base_proteins = read_base_proteins(data_folder/"base_proteins.csv")
settings = Settings.FromTOML(data_folder/"config.mainzer")
verbose = settings.settings["verbose"]

proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = run_lipido(
    mz=mz,
   intensity=intensity,
   base_proteins=base_proteins,
   base_lipids=base_lipids,
   **settings.settings)

analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')
output_folder = data_folder/spectrum_path.stem
final_folder = output_folder/analysis_time
final_folder.mkdir(parents=True, exist_ok=True)

proteins.to_csv(final_folder/"proteins.csv", index=False)
free_lipid_clusters.to_csv(final_folder/"free_lipid_clusters.csv", index=False)
centroids_df.to_csv(final_folder/"centroids.csv")
settings.save_toml(final_folder/"config.mainzer")

# plot_fit(centroids_df)