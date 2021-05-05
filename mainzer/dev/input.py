%load_ext autoreload
%autoreload 2

import json
import pandas as pd
import numpy as np
import tqdm

with open("mainzer/dev/settings.json") as f:
    settings = json.load(f)

from mainzer.isotope_ops import IsotopicEnvelopes
from mainzer.centroids import Centroids
from mainzer.graph_ops import get_regression_bigraph
from mainzer.models import fit_DeconvolvedUnderfit
from mainzer.timer import Timer
from mainzer.read import read
from mainzer.lipido import get_lipido_ions


molecules = pd.read_csv(settings['path_molecules'])
ions = get_lipido_ions(molecules, **settings)


deconvolve=True
isotopic_coverage=.95
isotopic_bin_size=.1
neighbourhood_thr=1.1
underfitting_quantile=0.05
fitting_to_void_penalty=1.0
verbose=True
verbose_output=True
pd.set_option('display.max_columns', None)


# /home/matteo/Projects/och_Kallol/unlipid/data/test/molecules.csv
# /home/matteo/Projects/och_Kallol/unlipid/data/07232020_Resolution50000_64bit.mzML
from mainzer.signal_ops import cluster_spectrum
from mainzer.centroids import Centroids
from mainzer.spectrum_preprocessing import get_mz2dmz, mz2dmz_ppm
import matplotlib.pyplot as plt


mz, intensity = read(settings['path_spectrum'])
clustered_spectrum_df = cluster_spectrum(mz, intensity)
clustered_spectrum_df.eval("""dmz = right_mz - left_mz
                           cluster_size = right_idx - left_idx+1
                           """, inplace=True)

orbi_res = pd.read_csv("../resolution_measurements_50000.csv")
mz2dmz_anirudha, mz2res, coefs = get_mz2dmz(orbi_res.mz, orbi_res.resolution)

b,a = coefs
a/-b

plt.scatter(orbi_res.mz, orbi_res.resolution)
xx = np.linspace(orbi_res.mz.min(), orbi_res.mz.max(), 100)
plt.plot(xx, mz2res(xx))
plt.show()

mzmz = np.linspace(mz.min(), mz.max(), 1000)

plt.plot(mzmz, mz2dmz_anirudha(mzmz))
plt.plot(mzmz, mz2dmz_ppm(mzmz))
plt.show()




plt.scatter(clustered_spectrum_df.mz_apex,
            clustered_spectrum_df.dmz,
            c=clustered_spectrum_df.cluster_size)
plt.plot(mzmz, mz2dmz_ppm(mzmz))
plt.show()
resolutions

# need filters here: or not!
clustered_spectrum_df["mz_right"] = np.minimum(clustered_spectrum_df.right_mz, 
                                               clustered_spectrum_df.mz_apex + 
                                                 mz2dmz_ppm(clustered_spectrum_df.mz_apex))
clustered_spectrum_df["mz_left"] = np.minimum(clustered_spectrum_df.left_mz, 
                                              clustered_spectrum_df.mz_apex - mz2dmz_ppm(clustered_spectrum_df.mz_apex))

min_intensity_threshold = 100
clustered_spectrum_df[clustered_spectrum_df.I_sum >= min_intensity_threshold]

centroids = Centroids(clustered_spectrum_df)
