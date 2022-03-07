from .matchmaker import Matchmaker


def single_precursor_regression(
    ions,
    centroids,
    isotopic_calculator,
    neighbourhood_thr: float=1.1,
    min_neighbourhood_intensity: float=100.0,
    max_expected_ppm_distance: float=15.0,
    underfitting_quantile: float=0.0,
    min_total_fitted_probability: float=.8,
    min_max_intensity_threshold: float=100,
    min_charge_sequence_length: int=1,
    verbose: bool=False,
    **kwargs
):
    matchmaker = Matchmaker(
        ions=ions,
        centroids=centroids,
        isotopic_calculator=isotopic_calculator,
    )
    matchmaker.get_isotopic_summaries()
    matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
    matchmaker.assign_isotopologues_to_centroids(min_neighbourhood_intensity, verbose)
    matchmaker.mark_ions_far_from_apexes_of_centroids(max_expected_ppm_distance)
    matchmaker.estimate_max_ion_intensity(underfitting_quantile)
    matchmaker.summarize_ion_assignments()
    matchmaker.get_theoretical_intensities_where_no_signal_was_matched()
    matchmaker.get_total_intensity_errors_for_maximal_intensity_estimates()
    matchmaker.mark_ions_with_probability_mass_outside_centroids(min_total_fitted_probability)
    matchmaker.mark_ions_with_low_maximal_intensity(min_max_intensity_threshold)
    matchmaker.mark_ions_not_in_charge_cluster(min_charge_sequence_length)
    return matchmaker


def turn_single_precursor_regression_chimeric(
    matchmaker,
    fitting_to_void_penalty=1.0, 
    merge_zeros=True,
    normalize_X=False,
    chimeric_regression_fits_cnt=3,
    min_chimeric_intensity_threshold=100,
    verbose=True,
    **kwargs
):
    if chimeric_regression_fits_cnt >= 2:
        print(f"Running chimeric_regression for the first time.")
    else:
        print(f"Running chimeric_regression for the first and final time.")
    matchmaker.build_regression_bigraph()
    matchmaker.fit_multiple_ion_regressions(fitting_to_void_penalty, 
                                            merge_zeros,
                                            normalize_X,
                                            verbose)
    for i in range(chimeric_regression_fits_cnt-1):
        if verbose:
            print(f"Running chimeric_regression for the {i+2} time.")
        matchmaker.filter_out_ions_with_low_chimeric_estimates(min_chimeric_intensity_threshold)
        matchmaker.reset_state_to_before_chimeric_regression()
        matchmaker.build_regression_bigraph()
        matchmaker.fit_multiple_ion_regressions(fitting_to_void_penalty, 
                                                merge_zeros,
                                                normalize_X,
                                                verbose)
    matchmaker.get_chimeric_groups()
    return matchmaker


def chimeric_regression(
    ions,
    centroids,
    isotopic_calculator,
    neighbourhood_thr=1.1,
    min_neighbourhood_intensity=100,
    max_expected_ppm_distance_to_cluster_apex=15,
    max_ppm_distance_top_prob_mz_cluster_apex=15,
    underfitting_quantile=0.0,
    min_total_fitted_probability=.8,
    min_max_intensity_threshold=100,
    min_charge_sequence_length=1,
    fitting_to_void_penalty=1.0,
    merge_zeros=True,
    normalize_X=False,
    chimeric_regression_fits_cnt=3,
    min_chimeric_intensity_threshold=100,
    verbose=False,
    **kwargs
):
    """Full metal chimeric regression."""
    matchmaker = single_precursor_regression(
        ions,
        centroids,
        isotopic_calculator,
        neighbourhood_thr,
        min_neighbourhood_intensity,
        max_ppm_distance,
        ppm_distance_type,
        underfitting_quantile,
        min_total_fitted_probability,
        min_max_intensity_threshold,
        min_charge_sequence_length,
        verbose
    )
    matchmaker = turn_single_precursor_regression_chimeric(
        matchmaker,
        fitting_to_void_penalty, 
        merge_zeros,
        normalize_X,
        min_chimeric_intensity_threshold,
        chimeric_regression_fits_cnt,
        verbose
    )
    return matchmaker


