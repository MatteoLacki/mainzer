from .matchmaker import Matchmaker


def single_precursor_regression(ions,
                                centroids,
                                isotopic_calculator,
                                neighbourhood_thr=1.1,
                                min_neighbourhood_intensity=100,
                                underfitting_quantile=0.0,
                                min_total_fitted_probability=.8,
                                min_max_intensity_threshold=100,
                                min_charge_sequence_length=1,
                                verbose=False):
    matchmaker = Matchmaker(ions,
                            centroids,
                            isotopic_calculator)
    matchmaker.get_isotopic_summaries()
    matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
    matchmaker.assign_isotopologues_to_centroids(min_neighbourhood_intensity, verbose)
    matchmaker.estimate_max_ion_intensity(underfitting_quantile)
    matchmaker.summarize_ion_assignments()
    matchmaker.add_ion_assignments_summary_to_ions()
    matchmaker.get_theoretical_intensities_where_no_signal_was_matched()
    matchmaker.get_total_intensity_errors_for_maximal_intensity_estimates()
    matchmaker.filter_ions_that_do_not_fit_centroids(min_total_fitted_probability)
    matchmaker.filter_ions_with_low_maximal_intensity(min_max_intensity_threshold)
    matchmaker.charge_ions_that_are_not_in_chargestate_sequence(min_charge_sequence_length)
    return matchmaker


def turn_single_precursor_regression_chimeric(matchmaker,
                                              fitting_to_void_penalty=1.0, 
                                              merge_zeros=True,
                                              normalize_X=False,
                                              chimeric_regression_fits_cnt=3,
                                              min_chimeric_intensity_threshold=100,
                                              verbose=True):
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


def chimeric_regression(ions,
                        centroids,
                        isotopic_calculator,
                        neighbourhood_thr=1.1,
                        min_neighbourhood_intensity=100,
                        underfitting_quantile=0.0,
                        min_total_fitted_probability=.8,
                        min_max_intensity_threshold=100,
                        min_charge_sequence_length=1,
                        fitting_to_void_penalty=1.0,
                        merge_zeros=True,
                        normalize_X=False,
                        chimeric_regression_fits_cnt=3,
                        min_chimeric_intensity_threshold=100,
                        verbose=False):
    """Full metal chimeric regression."""
    matchmaker = Matchmaker(ions, centroids, isotopic_calculator)
    matchmaker.get_isotopic_summaries()
    matchmaker.get_neighbourhood_intensities(neighbourhood_thr)
    matchmaker.assign_isotopologues_to_centroids(min_neighbourhood_intensity, verbose)
    matchmaker.estimate_max_ion_intensity(underfitting_quantile)
    matchmaker.summarize_ion_assignments()
    matchmaker.add_ion_assignments_summary_to_ions()
    matchmaker.get_theoretical_intensities_where_no_signal_was_matched()
    matchmaker.get_total_intensity_errors_for_maximal_intensity_estimates()
    matchmaker.filter_ions_that_do_not_fit_centroids(min_total_fitted_probability)
    matchmaker.filter_ions_with_low_maximal_intensity(min_max_intensity_threshold)
    matchmaker.charge_ions_that_are_not_in_chargestate_sequence(min_charge_sequence_length)
    matchmaker = turn_single_precursor_regression_chimeric(matchmaker,
                                                           fitting_to_void_penalty, 
                                                           merge_zeros,
                                                           normalize_X,
                                                           min_chimeric_intensity_threshold,
                                                           chimeric_regression_fits_cnt,
                                                           verbose)
    return matchmaker


