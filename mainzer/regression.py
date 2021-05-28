import numpy as np
import pandas as pd
import tqdm

from .isotope_ops import IsotopicEnvelopes
from .graph_ops import get_regression_bigraph
from .models import fit_DeconvolvedUnderfit

DEBUG = True


# ions_df = protein_mers
# centroids = filtered_centroids
def single_molecule_regression(centroids,
                               ions_df,
                               isotopic_coverage=.95,
                               isotopic_bin_size=.1,
                               neighbourhood_thr=1.1,
                               underfitting_quantile=0.05,
                               verbose=False,
                               simple_output=True,
                               **kwds):

    assert all(colname in ions_df.columns for colname in ("name","formula","charge")), "'ions_df' should be a pandas.DataFrame with columns 'name', 'formula', and 'charge'."
    assert isotopic_coverage >= 0 and isotopic_coverage <= 1, "'isotopic_coverage' is a probability and so must be between 0 and 1."
    assert underfitting_quantile >= 0 and underfitting_quantile <= 1, "'underfitting_quantile' is a probability and so must be between 0 and 1."

    if verbose:
        print("Getting Isotopic Envelopes")

    isotopic_envelopes = IsotopicEnvelopes(ions_df.formula.unique(),
                                           isotopic_coverage,
                                           isotopic_bin_size)
    # self = isotopic_envelopes
    # del self
    # formulas = ions_df.formula
    # charges = ions_df.charge
    # this part is totally weird...
    ions_df = ions_df.merge(isotopic_envelopes.charged_envelopes_summary(ions_df.formula, ions_df.charge))
    
    ion_idx = ['formula','charge']
    ions_df = ions_df.set_index(ion_idx)

    if verbose:
        print(f"Getting intensity of signal centroids in a {neighbourhood_thr}-neighbourhood of theoretical peaks.")

    minmax_signals = centroids.interval_query(ions_df.min_isospec_mz - neighbourhood_thr,
                                              ions_df.max_isospec_mz + neighbourhood_thr)
    
    ions_df["neighbourhood_intensity"] = minmax_signals.I_sum.groupby(ion_idx).sum()
    ions_df.neighbourhood_intensity.fillna(0, inplace=True)# NaN = 0 intensity
    del minmax_signals

    if verbose:
        print("Assigning isotopic envelope peaks to real signal centroids.")

    full_envelopes = isotopic_envelopes.to_frame(ions_df[ions_df.neighbourhood_intensity > 0].index)
    peak_assignments = centroids.point_query(full_envelopes.isospec_mz)
    peak_assignments = pd.concat([full_envelopes.iloc[peak_assignments.index],
                                  peak_assignments], axis=1)
    del full_envelopes

    assert np.all(peak_assignments.isospec_mz == peak_assignments.query_mz), "Something is wrong while creating peak assignments: query m/z not the same as isospec_mz."
    peak_assignments.drop("query_mz", axis=1, inplace=True) # "query_mz" is the same as "isospec_mz"
    peak_assignments = peak_assignments[peak_assignments.I_sum > 0]
    # need to have a measure of error

    # signal-based centroiding of isospec peaks! so we do not really need that binning!!!
    # watch out for summing "I_sum" twice! "I_sum" must be part of a group here:
    peak_assignments_clustered = peak_assignments.groupby(ion_idx + ["cluster","left_mz","right_mz","mz_apex","I_sum"]).agg({"isospec_prob":"sum", "isospec_mz":['min','max']}).reset_index()
    peak_assignments_clustered.columns = [f"{a}_{b}" if a == "isospec_mz" else a for a, b in peak_assignments_clustered.columns] # how I hate pandas..

    if verbose:
        print("Getting maximal intensity estimates.")

    max_intensity = peak_assignments_clustered.set_index(ion_idx)
    max_intensity = max_intensity.I_sum / max_intensity.isospec_prob
    max_intensity = max_intensity.groupby(ion_idx).quantile(underfitting_quantile)
    ions_df["maximal_intensity"] = max_intensity
    ions_df.maximal_intensity.fillna(0, inplace=True)

    if verbose:
        print("Summarizing assignments.")

    peak_assignments_summary = peak_assignments_clustered.groupby(ion_idx).agg({"isospec_prob":"sum", "I_sum":"sum", "cluster":"nunique"}).rename(columns={"isospec_prob":"isospec_prob_with_signal", "I_sum":"proximity_intensity", "cluster":"touched_centroids"})
    peak_assignments_summary["isospec_final_coverage"] = ions_df.isospec_final_coverage
    peak_assignments_summary["isospec_prob_without_signal"] = \
        peak_assignments_summary.isospec_final_coverage - \
        peak_assignments_summary.isospec_prob_with_signal

    # adding peak assignment summary to ions_df
    ions_df = pd.merge(ions_df,
                       peak_assignments_summary[['isospec_prob_with_signal',
                                                 'isospec_prob_without_signal',
                                                 'proximity_intensity',
                                                 'touched_centroids']], 
                       left_index=True, 
                       right_index=True,
                       how='left')
    
    ions_df.isospec_prob_with_signal.fillna(0, inplace=True)
    ions_df.proximity_intensity.fillna(0, inplace=True)
    ions_df.touched_centroids.fillna(0, inplace=True)
    ions_df.touched_centroids = ions_df.touched_centroids.astype(int)
    ions_df.isospec_prob_without_signal = np.where(ions_df.isospec_prob_without_signal.isna(), ions_df.isospec_final_coverage, ions_df.isospec_prob_without_signal)

    ions_df = ions_df.reset_index()

    if simple_output:
        return ions_df
    else:
        return ions_df, peak_assignments_clustered, peak_assignments_summary 


def single_molecule_regression_old(centroids,
                               ions_df,
                               isotopic_coverage=.95,
                               isotopic_bin_size=.1,
                               neighbourhood_thr=1.1,
                               underfitting_quantile=0.05,
                               verbose=False,
                               **kwds):

    assert all(colname in ions_df.columns for colname in ("name","formula","charge")), "'ions_df' should be a pandas.DataFrame with columns 'name', 'formula', and 'charge'."
    assert isotopic_coverage >= 0 and isotopic_coverage <= 1, "'isotopic_coverage' is a probability and so must be between 0 and 1."
    assert underfitting_quantile >= 0 and underfitting_quantile <= 1, "'underfitting_quantile' is a probability and so must be between 0 and 1."

    if verbose:
        print("Getting Isotopic Envelopes")

    isotopic_envelopes = IsotopicEnvelopes(ions_df.formula.unique(),
                                           isotopic_coverage,
                                           isotopic_bin_size)
    if DEBUG:
        print("Ctored")

    ions_df = ions_df.merge(isotopic_envelopes.charged_envelopes_summary(ions_df.formula, ions_df.charge))
    
    if DEBUG:
        print("merged")
    
    ion_idx = ['formula','charge']
    ions_df = ions_df.set_index(ion_idx)

    if verbose:
        print(f"Getting intensity of signal centroids in a {neighbourhood_thr}-neighbourhood of theoretical peaks.")

    minmax_signals = centroids.interval_query(ions_df.min_isospec_mz - neighbourhood_thr,
                                              ions_df.max_isospec_mz + neighbourhood_thr)
    ions_df["neighbourhood_intensity"] = minmax_signals.I_sum.groupby(ion_idx).sum()
    ions_df.neighbourhood_intensity.fillna(0, inplace=True)# NaN = 0 intensity
    del minmax_signals


    if verbose:
        print("Assigning isotopic envelope peaks to real signal centroids.")

    full_envelopes = isotopic_envelopes.to_frame(ions_df[ions_df.neighbourhood_intensity > 0].index)
    peak_assignments = centroids.point_query(full_envelopes.isospec_mz)
    peak_assignments = pd.concat([full_envelopes.iloc[peak_assignments.index],
                                  peak_assignments], axis=1)
    del full_envelopes

    assert np.all(peak_assignments.isospec_mz == peak_assignments.query_mz), "Something is wrong while creating peak assignments: query m/z not the same as isospec_mz."
    peak_assignments.drop("query_mz", axis=1, inplace=True) # "query_mz" is the same as "isospec_mz"
    peak_assignments = peak_assignments[peak_assignments.I_sum > 0]
    # need to have a measure of error

    # signal-based centroiding of isospec peaks! so we do not really need that binning!!!
    # watch out for summing "I_sum" twice! "I_sum" must be part of a group here:
    peak_assignments_clustered = peak_assignments.groupby(ion_idx + ["cluster","left_mz","right_mz","mz_apex","I_sum"]).agg({"isospec_prob":"sum", "isospec_mz":['min','max']}).reset_index()
    peak_assignments_clustered.columns = [f"{a}_{b}" if a == "isospec_mz" else a for a, b in peak_assignments_clustered.columns] # how I hate pandas..


    if verbose:
        print("Getting maximal intensity estimates.")

    max_intensity = peak_assignments_clustered.set_index(ion_idx)
    max_intensity = max_intensity.I_sum / max_intensity.isospec_prob
    max_intensity = max_intensity.groupby(ion_idx).quantile(underfitting_quantile)
    ions_df["maximal_intensity"] = max_intensity
    ions_df.maximal_intensity.fillna(0, inplace=True)


    if verbose:
        print("Summarizing assignments.")

    peak_assignments_summary = peak_assignments_clustered.groupby(ion_idx).agg({"isospec_prob":"sum", "I_sum":"sum", "cluster":"nunique"}).rename(columns={"isospec_prob":"isospec_prob_with_signal", "I_sum":"proximity_intensity", "cluster":"touched_centroids"})
    peak_assignments_summary["isospec_final_coverage"] = ions_df.isospec_final_coverage
    peak_assignments_summary["isospec_prob_without_signal"] = \
        peak_assignments_summary.isospec_final_coverage - \
        peak_assignments_summary.isospec_prob_with_signal

    # adding peak assignment summary to ions_df
    ions_df = pd.merge(ions_df,
                    peak_assignments_summary[['isospec_prob_with_signal',
                                              'isospec_prob_without_signal',
                                              'proximity_intensity',
                                              'touched_centroids']], 
                    left_index=True, 
                    right_index=True,
                    how='left')
    ions_df.isospec_prob_with_signal.fillna(0, inplace=True)
    ions_df.proximity_intensity.fillna(0, inplace=True)
    ions_df.touched_centroids.fillna(0, inplace=True)
    ions_df.touched_centroids = ions_df.touched_centroids.astype(int)
    ions_df.isospec_prob_without_signal = np.where(ions_df.isospec_prob_without_signal.isna(), ions_df.isospec_final_coverage, ions_df.isospec_prob_without_signal)

    return ( ions_df, 
             centroids,
             peak_assignments_clustered,
             peak_assignments_summary )



def multiple_molecule_regression(ions_df,
                                 centroids,
                                 peak_assignments_clustered, 
                                 peak_assignments_summary,
                                 fitting_to_void_penalty=1.0,
                                 verbose=False):
    ion_idx = ['formula','charge']

    if verbose:
        print("Building deconvolution graph.")

    G = get_regression_bigraph(peak_assignments_clustered, peak_assignments_summary)
    convoluted_ions_df = list(G.iter_convoluted_ions_df()) # groups of ions_df competing for signal explanation


    if verbose:
        print("Fitting multiple molecule regressions_df.")
    iter_problems = G.iter_regression_problems(merge_zeros=True,
                                               normalize_X=False) # think about this normalization

    if verbose:
        iter_problems = tqdm.tqdm(iter_problems, total=len(convoluted_ions_df))
    
    models = [fit_DeconvolvedUnderfit(X,Y, lam=fitting_to_void_penalty)
              for X,Y in iter_problems]

    estimates = pd.concat([m.coef for m in models])
    estimates.name = 'estimate'
    estimates.index.names = ion_idx

    ions_df["deconvolved_intensity"] = estimates
    ions_df.deconvolved_intensity.fillna(0, inplace=True)

    # getting errors
    peak_assignments_clustered = peak_assignments_clustered.merge(estimates,
                                                                  left_on=ion_idx,
                                                                  right_index=True).rename(columns={"estimate":"alpha"})
    peak_assignments_clustered.alpha.fillna(0, inplace=True)
    peak_assignments_clustered["under_estimate"] = \
        peak_assignments_clustered.alpha * peak_assignments_clustered.isospec_prob
    
    centroids.df["under_estimate"] = peak_assignments_clustered.groupby("cluster").under_estimate.sum()
    centroids.df.under_estimate.fillna(0.0, inplace=True)
    centroids.df["under_estimate_remainder"] = centroids.df.I_sum - centroids.df.under_estimate

    # assert np.all(ions_df.maximal_intensity >= ions_df.deconvolved_intensity + eps), "Maximal estimates lower than deconvolved!"
    #TODO: check if the difference is small: eps

    if verbose:
        print("Multiple regression done!")
    return ions_df, centroids, peak_assignments_clustered, G, convoluted_ions_df, models
