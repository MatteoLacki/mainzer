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
                        verbose=False):
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
    matchmaker.build_regression_bigraph()
    matchmaker.get_chimeric_ions()
    matchmaker.fit_multiple_ion_regressions(fitting_to_void_penalty, 
                                            merge_zeros,
                                            normalize_X,
                                            verbose)
    return matchmaker


def turn_single_precursor_regression_chimeric(matchmaker,
                                              fitting_to_void_penalty=1.0, 
                                              merge_zeros=True,
                                              normalize_X=False,
                                              verbose=True):
    matchmaker.build_regression_bigraph()
    matchmaker.get_chimeric_ions()
    matchmaker.fit_multiple_ion_regressions(fitting_to_void_penalty, 
                                            merge_zeros,
                                            normalize_X,
                                            verbose)
    return matchmaker
