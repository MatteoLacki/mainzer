import pathlib
from datetime import datetime
import json
import pandas as pd

from .read import read
from .deconv import estimate_intensities
from .baseline import strip_baseline
from .ion_generators import get_lipido_ions


def lipido_main(settings):
    analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')
    molecules = pd.read_csv(settings['path_molecules'])
    mz, intensity = read(settings['path_spectrum'])
#    output_folder = pathlib.Path(input("Ouput folder: ")).expanduser()
    output_folder = pathlib.Path(settings['output_folder'])
    output_folder.mkdir(parents=True, exist_ok=True)
    (output_folder/analysis_time).mkdir(parents=True, exist_ok=True)
    verbose = settings.settings["verbose"]

#    mz, intensity = mz[intensity > settings["min_intensity"]], intensity[intensity > settings["min_intensity"]]

    if verbose:
        print()
        print("Running Lipido with:")
        settings.print_summary()
        print()
        print("It's business time!")

        print("Getting ions")
    ions = get_lipido_ions(molecules, **(settings.settings))

    if verbose:
        print("Baseline correction")

    mz, intensity = strip_baseline(mz, intensity, settings["min_intensity"])
    if verbose:
        print("Estimating intenisities")

    ions, centroids, timings  = estimate_intensities(mz, intensity, ions, verbose_output=False, **(settings.settings))

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

    final_folder = output_folder/analysis_time

    if verbose:
        print("Saving results")

    proteins.to_csv(final_folder/"proteins.csv")
    lipid_clusters.to_csv(final_folder/"lipid_clusters.csv")
    centroids.df.to_csv(final_folder/"centroids.csv")

    with open(final_folder/"timings.json", "w") as jsonfile:
        json.dump(timings, jsonfile, indent=4)

    # stop outputting jsons?
    with open(final_folder/"settings.json", "w") as jsonfile:
        json.dump(settings.settings, jsonfile, indent=4)
    settings.save_toml(final_folder/"config.mainzer")

    if verbose:
        print("Thank you for letting Lipido do its job!")

