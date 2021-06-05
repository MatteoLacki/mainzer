%load_ext autoreload
%autoreload 2
from mainzer.lipido import lipido_IO
from mainzer.settings import Settings
from mainzer.plot import plot_fit

settings = Settings.FromTOML("settings.mainzer")
settings.settings["output_folder"] = "testres"

proteins, free_lipid_clusters, centroids_df = lipido_IO(settings)
plot_fit(centroids_df)