import numpy as np
import pandas as pd

import IsoSpecPy as isospec


class IsotopicEnvelopes(dict):
    def __init__(self,
                 formulas,
                 coverage=.95,
                 bin_size=.1,
                 get_minimal_pset=True,
                 PROTON_MASS=1.007276):
        self.coverage = coverage
        self.bin_size = bin_size
        self.get_minimal_pset = get_minimal_pset
        self.PROTON_MASS = PROTON_MASS
        for formula in formulas:
            iso = isospec.IsoTotalProb(prob_to_cover=self.coverage,
                                       formula=formula,
                                       get_minimal_pset=self.get_minimal_pset)
            if self.bin_size is not None:
                iso = iso.binned(self.bin_size)
            self[formula] = iso

    def formula2masses(self, formula):
        return self[formula].np_masses()

    def formula2mzs(self, formula, charge):
        return self[formula].np_masses() / charge + self.PROTON_MASS

    def formula2probs(self, formula):
        return self[formula].np_probs()

    def iter_envelope_summaries(self, formulas=None):
        for formula in formulas if formulas is not None else self:
            envelope = self[formula]
            masses = envelope.np_masses()
            probs = envelope.np_probs()
            yield {"formula": formula,
                   "isospec_peaks_count": len(envelope),
                   "isospec_final_coverage": probs.sum(),
                   "min_mass": masses.min(),
                   "max_mass": masses.max()}

    def envelopes_summary(self, formulas=None):
        return pd.DataFrame(self.iter_envelope_summaries(formulas))

    def charged_envelopes_summary(self, formulas, charges):
        res = self.envelopes_summary(formulas)
        res['charge'] = charges
        res.eval("""min_isospec_mz = min_mass / charge + @self.PROTON_MASS
                    max_isospec_mz = max_mass / charge + @self.PROTON_MASS""", inplace=True)
        res.drop(['min_mass','max_mass'], axis=1, inplace=True)
        return res

    def sizes(self):
        return {formula: len(iso) for formula,iso in self.items()}

    def to_frame(self, ions, return_arrays=False):
        theory_peaks_cnts = np.array([len(self[formula]) for formula, _ in ions])
        theory_peaks_cnt = theory_peaks_cnts.sum()
        mzs = np.zeros(theory_peaks_cnt, float)
        probs = np.zeros(theory_peaks_cnt, float)
        i_prev = 0
        for formula, charge in ions:
            iso = self[formula]
            i = i_prev + len(iso)
            mzs[i_prev:i] = iso.np_masses() / charge + self.PROTON_MASS
            probs[i_prev:i] = iso.np_probs()
            i_prev = i
        if return_arrays:
            return mzs, probs
        else:
            res = pd.DataFrame({"formula":np.repeat(ions.get_level_values("formula"), theory_peaks_cnts),
                                "charge": np.repeat(ions.get_level_values("charge"), theory_peaks_cnts),
                                "isospec_mz":mzs,
                                "isospec_prob":probs})
            return res


# def get_mzs_and_probs(formulas, charges, peak_counts, isotopic_envelopes, PROTON_MASS = 1.007276):
#     theory_peaks_cnt = peak_counts.sum()
#     mzs = np.zeros(theory_peaks_cnt, float)
#     probs = np.zeros(theory_peaks_cnt, float)
#     i_prev = 0
#     for formula, charge in zip(formulas, charges):
#         iso = isotopic_envelopes[formula]
#         i = i_prev + len(iso)
#         mzs[i_prev:i] = iso.np_masses() / charge + PROTON_MASS
#         probs[i_prev:i] = iso.np_probs()
#         i_prev = i
#     return mzs, probs


# def assign_peaks(formulas, charges, peak_counts, isotopic_envelopes, spectrum_intervals):
#     mzs, probs = get_mzs_and_probs(formulas, charges, peak_counts, isotopic_envelopes)

#     # not using values results in spread of the index from formulas... and a bug would follow!
#     peak_assignments = pd.DataFrame({"formula":np.repeat(formulas.values, peak_counts),
#                                      "charge": np.repeat(charges.values, peak_counts),
#                                      "mz":     mzs,
#                                      "prob":   probs })

#     theory_idxs, spectrum_intervals_idx = spectrum_intervals[peak_assignments.mz]
#     # making a series indexed with theory_idxs makes it easy to "merge by assignement" ..
#     spectrum_intervals_idx = pd.Series(spectrum_intervals_idx, index=theory_idxs)
#     # here! so that each intensity is properly assigned to its ion:
#     peak_assignments["spectrum_interval"] = spectrum_intervals_idx
#     peak_assignments.spectrum_interval = peak_assignments.spectrum_interval.fillna(-1).astype(int)
#     peak_assignments.eval("mzprob = mz*prob", inplace=True)
#     return peak_assignments