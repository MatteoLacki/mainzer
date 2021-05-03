import pathlib
from datetime import datetime
import json
import pandas as pd
from aa2atom import aa2atom, atom2str
from mainzer.read import read
from mainzer.deconv import estimate_intensities
from mainzer.settings import Settings
from mainzer.baseline import strip_baseline

aa2str = lambda aa: atom2str(aa2atom(aa))

from .molecule_ops import iter_mered_molecules, crosslink, molecules2df


def get_lipido_ions(molecules,
                    max_lipid_mers,
                    min_lipid_charge,
                    max_lipid_charge,
                    min_protein_charge,
                    max_protein_charge,
                    **kwargs):
    """Get a set of ions for lipido.

    Arguments:
        molecules (pd.DataFrame): A data frame with columns "group", "name", and "sequence_or_formula".
        max_lipid_mers (int): 
        min_lipid_charge (int):
        max_lipid_charge (int):
        min_protein_charge (int):
        max_protein_charge (int):
        **kwargs: other key-world arguments.
    Returns:
        pd.DataFrame: A data.frame with columns "name", "formula", and "charge": the general input for intensity estimation algorithms.
    """
    assert all(col in molecules.columns for col in ("group","name","sequence_or_formula")), "The csv with molecules should have columns 'group', 'name', and 'sequence'."
    proteins = molecules[molecules.group == "protein"]
    protein_formulas = []
    for seq in proteins.sequence_or_formula:
        formula = aa2str(seq)
        protein_formulas.append(formula)
        # except Exception:
        #     protein_formulas.append(seq)

    proteins = dict(zip(proteins.name, protein_formulas))

    lipids = molecules[molecules.group == "lipid"]
    lipids = dict(zip(lipids.name, lipids.sequence_or_formula))

    lipid_mers = dict(iter_mered_molecules(lipids, max_lipid_mers))
    lipid_protein_mers = dict(crosslink(lipid_mers, proteins))

    ions = pd.concat([molecules2df(lipid_mers, range(min_lipid_charge, max_lipid_charge+1)),
                      molecules2df(proteins, range(min_protein_charge, max_protein_charge+1)),
                      molecules2df(lipid_protein_mers, range(min_protein_charge, max_protein_charge+1))],
                      ignore_index=True)
    return ions


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
    with open(final_folder/"settings.json", "w") as jsonfile:
        json.dump(settings.settings, jsonfile, indent=4)
    settings.save_toml(final_folder/"config.mainzer")

    if verbose:
        print("Thank you for letting Lipido do its job!")

