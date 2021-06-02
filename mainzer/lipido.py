from datetime import datetime
import json
import pandas as pd
import pathlib

# from .baseline import strip_baseline
from .signal_ops import cluster_spectrum
from .isotope_ops import IsotopicCalculator
from .molecule_ops import mered_proteins, mered_lipids, molecules2df, crosslink, merge_free_lipids_and_promissing_proteins
from .read import read_spectrum, read_base_lipids, read_base_proteins
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

    proteins, free_lipid_clusters, centroids_df = \
        run_lipido(mz=mz,
                   intensity=intensity,
                   base_proteins=base_proteins,
                   base_lipids=base_lipids,
                   **settings.settings)

    final_folder = output_folder/analysis_time

    if verbose:
        print("Saving results.")
    proteins.to_csv(final_folder/"proteins.csv")
    free_lipid_clusters.to_csv(final_folder/"free_lipid_clusters.csv")
    centroids_df.to_csv(final_folder/"centroids.csv")
    settings.save_toml(final_folder/"config.mainzer")

    if verbose:
        print("Thank you for running Lipido!")

# defaults are set in settings.py! no need to repeat them here.
def run_lipido(# spectrum preprocessing
               mz,
               intensity,
               min_highest_intensity,
               min_mz,
               max_mz,
               # proteins
               base_proteins,
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
               # ion filters
               min_neighbourhood_intensity,
               min_charge_sequence_length,
               min_total_fitted_probability,
               # single molecule regression
               isotopic_coverage,
               isotopic_bin_size,
               neighbourhood_thr,
               underfitting_quantile,
               min_max_intensity_threshold,
               # multiple molecule regression
               run_chimeric_regression,
               fitting_to_void_penalty,
               # verbosity might be removed in favor of a logger
               verbose=False,
               **kwds):
    
    if verbose:
        print("Centroiding")#TODO: change for logger

    #TODO: reintroduce1 Michals baseline correction idea
    centroids_df = cluster_spectrum(mz, intensity)
    
    # this is simple: filerting on `max_intensity` in `centroids_df`
    filtered_centroids_df = centroids_df[(centroids_df.highest_intensity >= min_highest_intensity).values &\
                                         (centroids_df.left_mz >= min_mz).values &\
                                         (centroids_df.right_mz <= max_mz).values].copy()

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

    regression_kwds = {
        "centroids": filtered_centroids_df,
        "isotopic_calculator": isotopic_calculator,
        "neighbourhood_thr": neighbourhood_thr,
        "min_neighbourhood_intensity": min_neighbourhood_intensity,
        "underfitting_quantile": underfitting_quantile,
        "min_total_fitted_probability": min_total_fitted_probability,
        "min_max_intensity_threshold": min_max_intensity_threshold,
        "verbose": verbose
    }

    matchmaker = single_precursor_regression(
        protein_ions,
        min_charge_sequence_length=min_charge_sequence_length,
        **regression_kwds
    )
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
        print("Getting promissing protein-lipid centroids.")

    promissing_protein_lipid_complexes = \
        merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
                                                  promissing_protein_ions_df)
    
    protein_ions['contains_protein'] = True
    promissing_protein_lipid_complexes['contains_protein'] = True
    free_lipids_with_charge_df['contains_protein'] = False
    
    # Promissing Ions = Promissing Proteins + Protein Lipid Clusters + Free Lipid Clusters
    promissing_ions = pd.concat([protein_ions,
                                 promissing_protein_lipid_complexes,
                                 free_lipids_with_charge_df],
                                 ignore_index=True)

    full_matchmaker = single_precursor_regression(promissing_ions,
                                                  min_charge_sequence_length=1,
                                                  **regression_kwds)

    if run_chimeric_regression:
        if verbose:
            print("Performing chimeric regression.")
        full_matchmaker = turn_single_precursor_regression_chimeric(
            full_matchmaker,
            fitting_to_void_penalty, 
            merge_zeros=True,
            normalize_X=False,
            verbose=verbose
        )

    final_ions = full_matchmaker.ions.copy()
    column_order = ["name", "contains_protein", "formula", "charge"]
    sort_cols = ["charge","chimeric_intensity_estimate"] if run_chimeric_regression else ["charge","maximal_intensity_estimate"]
    final_ions = final_ions.sort_values(sort_cols, ascending=[True, False])    
    if run_chimeric_regression:
        column_order.append("chimeric_intensity_estimate")

    column_order += [ "maximal_intensity_estimate",
                      "neighbourhood_intensity",
                      "envelope_size",
                      "envelope_total_prob",
                      "envelope_min_mz",
                      "envelope_max_mz",
                      "envelope_proximity_intensity",
                      "envelope_matched_to_signal",
                      "envelope_unmatched_prob",
                      "explainable_centroids",
                      "single_precursor_fit_error",
                      "single_precursor_unmatched_estimated_intensity",
                      "single_precursor_total_error"]

    final_ions = final_ions[column_order]
    
    # make it easier to distinguish proteins..
    proteins = final_ions[final_ions.contains_protein].copy().drop(columns="contains_protein")
    free_lipid_clusters = final_ions[~final_ions.contains_protein].copy().drop(columns="contains_protein")
    
    centroids_df = pd.merge(
        centroids_df,
        full_matchmaker.centroids[["chimeric_intensity_in_centroid", "chimeric_remainder"]],
        left_index=True,
        right_index=True,
        how='left'
    )
    centroids_df.fillna(0, inplace=True)
    #TODO maybe it would make sense to 
    return proteins, free_lipid_clusters, centroids_df

