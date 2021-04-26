import json
import pandas as pd
import numpy as np

with open("mainzer/dev/settings.json") as f:
    settings = json.load(f)

from mainzer.isotope_ops import IsotopicEnvelopes
from mainzer.centroids import Centroids
from mainzer.graph_ops import get_regression_bigraph
from mainzer.models import fit_DeconvolvedUnderfit
from mainzer.timer import Timer
from mainzer.read import read
from mainzer.lipido import get_lipido_ions

mz, intensity = read(settings['path_spectrum'])
molecules = pd.read_csv(settings['path_molecules'])
ions = get_lipido_ions(molecules, **settings)


deconvolve=True
isotopic_coverage=.95
isotopic_bin_size=.1
neighbourhood_thr=1.1
underfitting_quantile=0.05
fitting_to_void_penalty=1.0
verbose=False
verbose_output=True

pd.set_option('display.max_columns', None)
