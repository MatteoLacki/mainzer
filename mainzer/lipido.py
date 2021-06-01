from datetime import datetime
import json
import pandas as pd
import pathlib

from .baseline import strip_baseline
from .signal_ops import cluster_spectrum
from .isotope_ops import IsotopicCalculator
from .regression import single_molecule_regression, multiple_molecule_regression
from .molecule_ops import mered_proteins, mered_lipids, molecules2df, crosslink, merge_free_lipids_and_promissing_proteins
from .read import read_spectrum, read_molecules_for_lipido
from .regression import single_precursor_regression, turn_single_precursor_regression_chimeric


def lipido_IO(settings):
    """ 
    Read in necessary data and saves the outputs.
    """
    analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')

    base_lipids   = read_base_lipids(settings['path_base_lipids'])
    base_proteins = read_base_proteins(settings['path_base_proteins'])
    mz, intensity = read_spectrum(settings['path_spectrum'])
    
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

    # proteins, lipid_clusters, centroids = run_lipido(mz=mz,
    #                                                  intensty=intensity,
    #                                                  base_proteins=base_proteins,
    #                                                  base_lipids=base_lipids,
    #                                                  **settings.settings)

    # final_folder = output_folder/analysis_time

    # if verbose:
    #     print("Saving results.")
    # proteins.to_csv(final_folder/"proteins.csv")
    # lipid_clusters.to_csv(final_folder/"lipid_clusters.csv")
    # centroids.df.to_csv(final_folder/"centroids.csv")
    # settings.save_toml(final_folder/"config.mainzer")

    # if verbose:
    #     print("Thank you for letting Lipido do its job!")


# defaults are set in settings.py! no need to repeat them here.
def run_lipido(# spectrum preprocessing
               mz,
               intensity,
               min_intensity_threshold,
               min_mz,
               max_mz,
               # proteins
               base_proteins,
               min_protein_mers,
               max_protein_mers,
               min_protein_cluster_charge,
               max_protein_cluster_charge,
               # lipids
               base_lipids,
               min_lipid_mers,
               max_lipid_mers,
               min_free_lipid_cluster_charge,
               max_free_lipid_cluster_charge,
               # crosslinking
               only_heteromers, # discard homomers
               # charge based filter: turn off by setting to 1
               min_charge_sequence_length,
               # single molecule regression
               isotopic_coverage,
               isotopic_bin_size,
               neighbourhood_thr,
               underfitting_quantile,
               minimal_maximal_intensity_threshold,
               # multiple molecule regression
               deconvolve,
               fitting_to_void_penalty,
               # verbosity might be removed in favor of a logger
               verbose=False,
               **kwds):
    
    if verbose:
        print("Centroiding")#TODO: change for logger

    #TODO: reintroduce Michals baseline correction idea
    cenroids_df = cluster_spectrum(mz, intensity)
    
    # this is simple: filerting on `max_intensity` in `cenroids_df`
    filtered_cenroids_df = cenroids_df[(cenroids_df.I_max >= min_intensity_threshold).values &\
                                       (cenroids_df.left_mz >= min_mz).values &\
                                       (cenroids_df.right_mz <= max_mz).values].copy()

    if verbose:
        print("Getting proteins mers")#TODO: change for logger
    # initially we search for proteins only

    protein_mers = mered_proteins(base_proteins, only_heteromers)
    protein_ions = molecules2df(protein_mers,
                                range(min_protein_cluster_charge,
                                      max_protein_cluster_charge+1))

    if verbose:
        print("Setting up IsoSpec (soon to defeat NeutronStar, if it already does not).")
    isotopic_calculator = IsotopicCalculator(isotopic_coverage, isotopic_bin_size)

    if verbose:
        print("Checking for promissing proteins")

    matchmaker = single_precursor_regression(protein_ions,
                                             filtered_cenroids_df,
                                             isotopic_calculator,
                                             neighbourhood_thr,
                                             underfitting_quantile,
                                             min_total_fitted_probability,
                                             min_max_intensity_threshold,
                                             min_charge_sequence_length)
    promissing_protein_ions_df = matchmaker.ions
    
    if verbose:
        print("Getting lipid mers")
    
    free_lipid_mers = dict(mered_lipids(base_lipids,
                                        min_lipid_mers,
                                        max_lipid_mers))
    free_lipids_no_charge_df = molecules2df(free_lipid_mers)
    free_lipids_with_charge_df = \
        molecules2df(free_lipid_mers, 
                     charges = range(min_free_lipid_cluster_charge,
                                     max_free_lipid_cluster_charge))

    if verbose:
        print("Getting promissing protein-lipid cenroids.")

    promissing_protein_lipid_complexes = \
        merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
                                                  promissing_protein_ions_df)
    
    # Promissing Ions = Promissing Proteins + Protein Lipid Clusters + Free Lipid Clusters
    promissing_ions = pd.concat([protein_ions,
                                 promissing_protein_lipid_complexes,
                                 free_lipids_with_charge_df],
                                 ignore_index=True)

    full_matchmaker = single_precursor_regression(promissing_ions,
                                                  filtered_cenroids_df,
                                                  isotopic_calculator,
                                                  neighbourhood_thr,
                                                  underfitting_quantile,
                                                  min_total_fitted_probability,
                                                  min_max_intensity_threshold,
                                                  min_charge_sequence_length)

    if deconvolve:
        if verbose:
            print("Performing chimeric regression.")
        full_matchmaker = turn_single_precursor_regression_chimeric(
            full_matchmaker,
            fitting_to_void_penalty, 
            merge_zeros=True,
            normalize_X=False,
            verbose=verbose
        )

    #TODO: ??? additional tagging of ions based on found results.?????



    column_order = ["name"]
    if "deconvolved_intensity" in ions.columns:
        ions = ions.sort_values(['charge','deconvolved_intensity'], ascending=[True, False])
        column_order.append("deconvolved_intensity")
    else:
        ions = ions.sort_values(['charge','maximal_intensity'], ascending=[True, False])

    column_order.extend([ "maximal_intensity",
                          "proximity_intensity",
                          "neighbourhood_intensity",
                          "isospec_final_coverage",
                          "isospec_prob_with_signal",
                          "isospec_prob_without_signal",
                          "touched_centroids",
                          "isospec_peaks_count",
                          "min_isospec_mz",
                          "max_isospec_mz" ])

    ions = ions[column_order]
    protein_names = "|".join(molecules[molecules.group == "protein"].name)
    (_, lipid_cenroids), (_, proteins) = ions.groupby(ions.name.str.contains(protein_names), group_keys=False)

    return proteins, lipid_cenroids, centroids


