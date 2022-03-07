import datetime
import math
import pandas as pd
import pathlib

# from .baseline import strip_baseline
import mainzer.data_frame_ops
import mainzer.isotope_ops
import mainzer.molecule_ops
import mainzer.read
import mainzer.regression
import mainzer.signal_ops


def lipido_IO(settings, output_folder):
    """ 
    Read in necessary data and saves the outputs.
    """
    analysis_time = datetime.datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')

    base_lipids   = mainzer.read.read_base_lipids(settings['path_base_lipids'])
    base_proteins = mainzer.read.read_base_proteins(settings['path_base_proteins'])
    print(settings['path_spectrum'])
    mz, intensity = mainzer.read.read_spectrum(settings['path_spectrum'])
    
    output_folder = pathlib.Path(output_folder)
    final_folder = output_folder/analysis_time
    
    verbose = settings.settings["verbose"]

    if verbose:
        print()
        print("Running Lipido with:")
        settings.print_summary()
        print()
        print("It's business time!")

    proteins, free_lipid_clusters, simple_proteins, simple_free_lipid_clusters, centroids_df = \
        run_lipido(mz=mz,
                   intensity=intensity,
                   base_proteins=base_proteins,
                   base_lipids=base_lipids,
                   params=settings.settings,
                   verbose=verbose)

    if verbose:
        print("Saving results.")

    final_folder.mkdir(parents=True, exist_ok=True)

    simple_proteins.to_csv(final_folder/"simple_protein_report.csv",
                           index=False)
    simple_free_lipid_clusters.to_csv(final_folder/"simple_free_lipid_clusters_report.csv",
                                      index=False)
    proteins.to_csv(final_folder/"extended_proteins_report.csv",
                    index=False)
    free_lipid_clusters.to_csv(final_folder/"extended_free_lipid_clusters_report.csv",
                               index=False)
    centroids_df.to_csv(final_folder/"centroids.csv")
    settings.save_toml(final_folder/"config.mainzer")

    if verbose:
        print("Thank you for running Lipido!")

    return proteins, free_lipid_clusters, centroids_df


def run_lipido(mz,
               intensity,
               base_proteins,
               base_lipids,
               params,
               verbose=False,
               debug=False):
    if verbose:
        print("Centroiding")#TODO: change for a logger... wait, whaaat? No way we will ever get so orderly.

    #TODO: reintroduce Michals baseline correction idea: the intensity correction might be dependent upon
    # mass to charge ratio.
    centroids_df = mainzer.signal_ops.cluster_spectrum(mz, intensity)
    
    # this is simple: filerting on `max_intensity` in `centroids_df`
    # do not use df.query: does not work on Windows after py2exe
    filtered_centroids_df = centroids_df[
        (centroids_df.highest_intensity >= params["min_highest_intensity"]).values &\
        (centroids_df.left_mz           >= params["min_mz"]).values &\
        (centroids_df.right_mz          <= params["max_mz"]).values
    ].copy()

    if verbose:
        print("Getting proteins mers")#TODO: change for logger
    # initially we search for proteins only

    protein_mers = mainzer.molecule_ops.mered_proteins(base_proteins, params["only_heteromers"])
    protein_ions = mainzer.molecule_ops.molecules2df(
        molecules=protein_mers,
        charges=range(params["min_protein_cluster_charge"], params["max_protein_cluster_charge"]+1)
    )

    if verbose:
        print("Setting up IsoSpec (soon to defeat NeutronStar, if it already does not).")
    isotopic_calculator = mainzer.isotope_ops.IsotopicCalculator(
        coverage=params["isotopic_coverage"], 
        bin_size=params["isotopic_bin_size"]
    )

    if verbose:
        print("Checking for promissing proteins")

    regression_kwds = params.copy()
    regression_kwds["centroids"] = filtered_centroids_df
    regression_kwds["isotopic_calculator"] = isotopic_calculator
    regression_kwds["verbose"] = verbose
    del regression_kwds["min_charge_sequence_length"]

    protein_ions_matchmaker = mainzer.regression.single_precursor_regression(
        ions=protein_ions,
        min_charge_sequence_length=params["min_charge_sequence_length"],
        **regression_kwds
    )
    promissing_protein_ions_df = protein_ions_matchmaker.ions[protein_ions_matchmaker.ions.reason_for_filtering_out=="none"].copy()
    
    if verbose:
        print("Getting lipid mers")
    
    free_lipid_mers = dict(
        mainzer.molecule_ops.mered_lipids(
            base_lipids=base_lipids,
            min_lipid_mers=params["min_lipid_mers"],
            max_lipid_mers=params["max_lipid_mers"]
        )
    )
    free_lipids_no_charge_df = mainzer.molecule_ops.molecules2df(free_lipid_mers)
    free_lipids_with_charge_df = mainzer.molecule_ops.molecules2df(
        free_lipid_mers, 
        charges=range(params["min_free_lipid_cluster_charge"],
                      params["max_free_lipid_cluster_charge"])
    )

    if verbose:
        print("Getting promissing protein-lipid centroids.")

    promissing_protein_lipid_complexes = mainzer.molecule_ops.merge_free_lipids_and_promissing_proteins(
        free_lipids_no_charge_df,
        promissing_protein_ions_df,
    )
    
    promissing_proteins = promissing_protein_ions_df[["name","formula","charge"]].copy()
    promissing_proteins['contains_protein'] = True
    promissing_protein_lipid_complexes['contains_protein'] = True
    free_lipids_with_charge_df['contains_protein'] = False
    
    # Promissing Ions = Promissing Proteins + Protein Lipid Clusters + Free Lipid Clusters
    promissing_ions = pd.concat(
        [
            promissing_proteins,
            promissing_protein_lipid_complexes,
            free_lipids_with_charge_df
        ],
        ignore_index=True
    )

    full_matchmaker = mainzer.regression.single_precursor_regression(
        promissing_ions,
        min_charge_sequence_length=1,
        **regression_kwds
    )

    run_chimeric_regression = params["chimeric_regression_fits_cnt"] > 0

    if run_chimeric_regression:
        if verbose:
            print("Performing chimeric regression.")
        full_matchmaker = mainzer.regression.turn_single_precursor_regression_chimeric(
            full_matchmaker,
            params["fitting_to_void_penalty"], 
            merge_zeros=True,
            normalize_X=False,#TODO: what about this??? Should it or not?
            chimeric_regression_fits_cnt=params["chimeric_regression_fits_cnt"],
            min_chimeric_intensity_threshold=params["min_chimeric_intensity_threshold"],
            verbose=verbose
        )

    #TODO put all this into Matchmaker 
    final_ions = full_matchmaker.ions[full_matchmaker.ions.reason_for_filtering_out == "none"].copy()

    # TODO: rething the settings... the column names should be there
    if params["rounding"] != -1:
        mz_round = math.ceil(-math.log10(params["isotopic_bin_size"]))
        columns2int = [
            "neighbourhood_intensity",
            "maximal_intensity_estimate",
            "single_precursor_fit_error",
            "envelope_proximity_intensity",
            "single_precursor_unmatched_estimated_intensity",
            "single_precursor_total_error"
        ]
        if run_chimeric_regression:
            columns2int.append("chimeric_intensity_estimate")
        final_ions = mainzer.data_frame_ops.round_df(
            final_ions,
            {"envelope_total_prob":          params["rounding"],
             "envelope_min_mz":              mz_round,
             "envelope_max_mz":              mz_round,
             "envelope_top_prob_mz":         mz_round,
             "envelope_matched_to_signal":   params["rounding"],
             "envelope_unmatched_prob":      params["rounding"]},
            columns2int
        )

    column_order = ["name", "contains_protein", "formula", "charge"]
    sort_cols = ["charge","chimeric_intensity_estimate"] if run_chimeric_regression else ["charge","maximal_intensity_estimate"]
    final_ions = final_ions.sort_values(sort_cols, ascending=[True, False])    
    if run_chimeric_regression:
        column_order.append("chimeric_intensity_estimate")

    column_order += [ "maximal_intensity_estimate",
                      "neighbourhood_intensity",
                      "maximal_intensity_estimate",
                      "neighbourhood_intensity"]
    if run_chimeric_regression:
        column_order.append("chimeric_group")
    column_order += [ "envelope_size",
                      "envelope_total_prob",
                      "envelope_min_mz",
                      "envelope_max_mz",
                      "envelope_top_prob_mz",
                      "envelope_proximity_intensity",
                      "envelope_matched_to_signal",
                      "envelope_unmatched_prob",
                      "explainable_centroids",
                      "single_precursor_fit_error",
                      "single_precursor_unmatched_estimated_intensity",
                      "single_precursor_total_error"]

    final_ions = final_ions[column_order]

    simple_ions_cols = [
        "name",
        "charge",
        "envelope_top_prob_mz",
        "chimeric_intensity_estimate" if run_chimeric_regression else "maximal_intensity_estimate"
    ]
    if run_chimeric_regression:
        simple_ions_cols.append("chimeric_group")

    simple_ions = final_ions[simple_ions_cols]
    simple_ions = simple_ions.rename(columns={"envelope_top_prob_mz":"top_probable_mz"})

    simple_proteins = simple_ions[final_ions.contains_protein]
    simple_free_lipid_clusters = simple_ions[~final_ions.contains_protein]

    # make it easier to distinguish proteins..
    proteins = final_ions[final_ions.contains_protein].copy().drop(columns="contains_protein")
    free_lipid_clusters = final_ions[~final_ions.contains_protein].copy().drop(columns="contains_protein")
    
    if run_chimeric_regression:
        centroids_df = pd.merge(
            centroids_df,
            full_matchmaker.centroids[["chimeric_intensity_in_centroid", "chimeric_remainder"]],
            left_index=True,
            right_index=True,
            how='left'
        )
        centroids_df.fillna(0, inplace=True)

        if params["rounding"] != -1:
            centroids_df = mainzer.data_frame_ops.round_df(
                centroids_df,
                {"mz_apex":     mz_round,
                 "left_mz":     mz_round,
                 "right_mz":    mz_round},
                ["highest_intensity",
                 "integrated_intensity",
                 "chimeric_intensity_in_centroid",
                 "chimeric_remainder"]
            )
        centroids_df.drop(columns=["left_idx", "right_idx"], inplace=True)
    
    if debug:
        return (
            proteins,
            free_lipid_clusters,
            simple_proteins,
            simple_free_lipid_clusters,
            centroids_df,
            full_matchmaker,
            protein_ions_matchmaker
        )
    else:
        return (
            proteins,
            free_lipid_clusters,
            simple_proteins,
            simple_free_lipid_clusters,
            centroids_df
        )

