%load_ext autoreload
%autoreload 2
import pathlib
import pandas as pd
pd.options.display.max_columns = None
import pprint

import mainzer.read
import mainzer.plot
import mainzer.settings
import mainzer.lipido
import mainzer.data_frame_ops
import mainzer.isotope_ops
import mainzer.molecule_ops
import mainzer.regression
import mainzer.signal_ops



data_folder = pathlib.Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/small_spectrum")
# data_folder = pathlib.Path("/home/matteo/Projects/och_Kallol/mainzer/test_data/bigger_spectrum")

if data_folder.exists():
    print(f"Folder {data_folder} exists.")

spectrum_path = next(data_folder.glob("*.mzML"))
mz, intensity = mainzer.read.read_spectrum(str(spectrum_path))
# mainzer.plot.plot_spectrum(mz, intensity)

base_lipids = mainzer.read.read_base_lipids(data_folder/"base_lipids.csv")
base_proteins = mainzer.read.read_base_proteins(data_folder/"base_proteins.csv")
settings = mainzer.settings.Settings.FromTOML(data_folder/"config.mainzer")


params = settings.settings
pprint.pprint(params)
verbose=True

proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = mainzer.lipido.run_lipido(
    mz=mz,
   intensity=intensity,
   base_proteins=base_proteins,
   base_lipids=base_lipids,
   params=settings.settings,
   verbose=True)

mainzer.plot.plot_fit(centroids_df)
free_lipid_clusters.head()
simple_proteins.head()
simple_free_lipid_clusters
centroids_df.head()
