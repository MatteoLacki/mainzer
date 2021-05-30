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
from .molecule_filters import charge_sequence_filter
from .intervals import IntervalQuery, IonsCentroids


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
        print("Centroiding") #TODO: change for logger

    #TODO: reintroduce Michals baseline correction idea
    clusters_df = cluster_spectrum(mz, intensity)
    
    # this is simple: filerting on `max_intensity` in `clusters_df`
    filtered_clusters_df = clusters_df[(clusters_df.I_max >= min_intensity_threshold).values &\
                                       (clusters_df.left_mz >= min_mz).values &\
                                       (clusters_df.right_mz <= max_mz).values].copy()

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

    ions_centroids = IonsCentroids(protein_ions,
                                   filtered_clusters_df,
                                   isotopic_calculator)
    ions_centroids.get_isotopic_summaries()
    


    # ions_df = protein_ions.copy()
    # centroids = filtered_clusters_df
    protein_ions = single_molecule_regression(filtered_clusters_df,
                                              protein_ions,
                                              isotopic_calculator,
                                              neighbourhood_thr,
                                              underfitting_quantile,
                                              verbose)



    # `protein_ions` is updated by regression results
    
    promissing_protein_ions_df = \
        protein_ions[protein_ions.maximal_intensity >= minimal_maximal_intensity_threshold].copy()

    #TODO: ERROR! ERROR! The esimates of maximal intensity should not exceed the spectral intensity. Something is massively wrong now.
    #TODO: change the "simple_molecule_regression" into a class that could output errors and report on different metrics.
    import matplotlib.pyplot as plt
    plt.scatter(
        promissing_protein_ions_df.isospec_prob_with_signal,
        np.log10(promissing_protein_ions_df.maximal_intensity/promissing_protein_ions_df.neighbourhood_intensity))
    plt.show()
    promissing_protein_ions_df.query("neighbourhood_intensity<maximal_intensity")

    
    if min_charge_sequence_length > 1:
        promissing_protein_ions_df = charge_sequence_filter(promissing_protein_ions_df,
                                                            min_charge_sequence_length)
    #TODO: might introduce extra criteria to filter more things out

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
        print("Getting promissing protein-lipid clusters.")

    promissing_protein_lipid_complexes = \
        merge_free_lipids_and_promissing_proteins(free_lipids_no_charge_df,
                                                  promissing_protein_ions_df)
    
    # (promissing) Protein Lipid Clusters + Free Lipid Clusters
    PLC_FLC = pd.concat([promissing_protein_lipid_complexes,
                         free_lipids_with_charge_df],
                         ignore_index=True)

    # ions_df need to have fresh integer key... I am sorry Michał, I created a monster.
    PLC_FLC = single_molecule_regression(filtered_centroids,
                                         PLC_FLC,
                                         isotopic_coverage,
                                         isotopic_bin_size,
                                         neighbourhood_thr,
                                         underfitting_quantile,
                                         verbose)
    
    PLC_FLC_intense = PLC_FLC[ PLC_FLC.maximal_intensity >= minimal_maximal_intensity_threshold].copy()

    import matplotlib.pyplot as plt
    
    plt.hist(PLC_FLC_intense.isospec_prob_with_signal)
    plt.show()
    plt.scatter(np.log(PLC_FLC_intense.maximal_intensity),
                np.log(PLC_FLC_intense.proximity_intensity),
                c=PLC_FLC_intense.isospec_prob_with_signal,
                s=1)
    plt.show()
    plt.scatter(
                PLC_FLC_intense.isospec_prob_with_signal,
                np.log10(PLC_FLC_intense.maximal_intensity / PLC_FLC_intense.neighbourhood_intensity),
                s=1)
    plt.show()
    # I need to export some misfits information.

    # run full


    # dict(crosslink(protein_mers, free_lipid_mers)) # too big!

    # if verbose:
    #     print("Looking for charged proteins.")
    
    
    

    # ions, centroids, peak_assignments_clustered, peak_assignments_summary = \
    #     single_molecule_regression( centroids,
    #                                 ions,
    #                                 isotopic_coverage,
    #                                 isotopic_bin_size,
    #                                 neighbourhood_thr,
    #                                 underfitting_quantile,
    #                                 verbose )
    

    # if verbose:
    #     print("Performing singe molecule deconvolution.")
    # ions, centroids, peak_assignments_clustered, peak_assignments_summary = \
    #     single_molecule_regression( centroids,
    #                                 ions,
    #                                 isotopic_coverage,
    #                                 isotopic_bin_size,
    #                                 neighbourhood_thr,
    #                                 underfitting_quantile,
    #                                 verbose )
    
    # here we need some filtering of ions

    if deconvolve:
        if verbose:
            print("Performing multiple molecule deconvolution.")
        ions, centroids, peak_assignments_clustered, G, convoluted_ions, models = \
            multiple_molecule_regression( ions,
                                          centroids,
                                          peak_assignments_clustered, 
                                          peak_assignments_summary,
                                          fitting_to_void_penalty,
                                          verbose )

    #TODO: additional tagging of ions based on found results.

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
    (_, lipid_clusters), (_, proteins) = ions.groupby(ions.name.str.contains(protein_names), group_keys=False)

    return proteins, lipid_clusters, centroids


