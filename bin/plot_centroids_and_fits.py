import argparse
import pandas as pd
import pathlib

from mainzer.plot import plot_fit

ap = argparse.ArgumentParser(description='Plot centroided spectrum with fits.',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('path_to_centroids_csv',
                help="Path to configuration file.",
                type=pathlib.Path)
ap = ap.parse_args()

centroids_df = pd.read_csv(ap.path_to_centroids_csv)
plot_fit(centroids_df)