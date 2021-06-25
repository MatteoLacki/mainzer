import pandas as pd
import pathlib
import sys

from .lipido import run_lipido
from .read import \
    read_spectrum,\
    read_base_lipids,\
    read_base_proteins
from .settings import Settings

test_data = pathlib.Path("test_data")

def test_gua():
    spectrum_path = test_data/"small_spectrum"/"09142020_SV_VAMP2_after_ON_Dialysis_2ndRoundCentrifugation-qb.mzML"

    mz, intensity = read_spectrum(spectrum_path)
    settings = Settings.FromTOML(spectrum_path.parent/"config.mainzer")
    base_lipids   = read_base_lipids(settings['path_base_lipids'])
    base_proteins = read_base_proteins(settings['path_base_proteins'])

    proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = \
        run_lipido(mz=mz,
                   intensity=intensity,
                   base_proteins=base_proteins,
                   base_lipids=base_lipids,
                   **settings.settings)

    old_proteins = pd.read_csv("")
    old_free_lipid_clusters = pd.read_csv("")
    old_settings = Settings.FromTOML()