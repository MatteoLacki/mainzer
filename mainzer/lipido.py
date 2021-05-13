from datetime import datetime
import json
import pandas as pd
import pathlib

from .baseline import strip_baseline
from .deconv import estimate_intensities
from .ion_generators import get_lipido_ions
from .read import read
from .centroiding import centroid_spectrum


def lipido_IO(settings):
    """ 
    Read in necessary data and saves the outputs.
    """
    analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')
    molecules = pd.read_csv(settings['path_molecules'])
    mz, intensity = read(settings['path_spectrum'])
    output_folder = pathlib.Path(settings['output_folder'])
    output_folder.mkdir(parents=True, exist_ok=True)
    (output_folder/analysis_time).mkdir(parents=True, exist_ok=True)
    verbose = settings.settings["verbose"]

    if verbose:
        print()
        print("Running Lipido with:")
        settings.print_summary()
        print()
        print("It's business time!")

    proteins, lipid_clusters, centroids = \
        run_lipido(mz, intensity, molecules, verbose,
                   settings_dict=settings.settings)

    final_folder = output_folder/analysis_time

    if verbose:
        print("Saving results.")
    proteins.to_csv(final_folder/"proteins.csv")
    lipid_clusters.to_csv(final_folder/"lipid_clusters.csv")
    centroids.df.to_csv(final_folder/"centroids.csv")
    settings.save_toml(final_folder/"config.mainzer")

    if verbose:
        print("Thank you for letting Lipido do its job!")



def run_lipido(mz, intensity, molecules, verbose, settings_dict):
    #TODO: wrap functions into some logger to get this right and time it.
    if verbose:
        print("Getting ions")
    ions = get_lipido_ions(molecules, **settings_dict)

    if verbose:
        print("Baseline correction")
    # mz, intensity = strip_baseline(mz, intensity, settings["min_intensity"])    

    if verbose:
        print("Centroiding")
    centroids = centroid_spectrum(mz, intensity)

    if verbose:
        print("Estimating intensities")
    ions, centroids = estimate_intensities(centroids,
                                           ions,
                                           **settings_dict)

    column_order = ["name"]
    if "deconvolved_intensity" in ions.columns:
        ions = ions.sort_values(['charge','deconvolved_intensity'], ascending=[True, False])
        column_order.append("deconvolved_intensity")
    else:
        ions = ions.sort_values(['charge','maximal_intensity'], ascending=[True, False])

    column_order.extend(["maximal_intensity",
                         "proximity_intensity",
                         "neighbourhood_intensity",
                         "isospec_final_coverage",
                         "isospec_prob_with_signal",
                         "isospec_prob_without_signal",
                         "touched_centroids",
                         "isospec_peaks_count",
                         "min_isospec_mz",
                         "max_isospec_mz"])

    ions = ions[column_order]
    protein_names = "|".join(molecules[molecules.group == "protein"].name)
    (_, lipid_clusters), (_, proteins) = ions.groupby(ions.name.str.contains(protein_names), group_keys=False)

    return proteins, lipid_clusters, centroids


