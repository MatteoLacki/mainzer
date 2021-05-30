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

    def iter_envelope_summaries(self):
        for formula in self:
            envelope = self[formula]
            masses = envelope.np_masses()
            probs = envelope.np_probs()
            yield {"formula": formula,
                   "isospec_peaks_count": len(envelope),
                   "isospec_final_coverage": probs.sum(),
                   "min_mass": masses.min(),
                   "max_mass": masses.max()}

    def envelopes_summary(self):
        return pd.DataFrame(self.iter_envelope_summaries())

    def ions_summary(self, ions_df, copy_ions_df=True):
        assert all(col in ions_df.columns for col in ("formula", "charge"))
        if copy_ions_df:
            ions_df = ions_df.copy()
        res = self.envelopes_summary()
        ions_df = ions_df.merge(res)
        ions_df["isospec_min_mz"] = ions_df.min_mass / ions_df.charge + self.PROTON_MASS
        ions_df["isospec_max_mz"] = ions_df.max_mass / ions_df.charge + self.PROTON_MASS
        ions_df.drop(['min_mass','max_mass'], axis=1, inplace=True)
        return ions_df

    def sizes(self):
        return {formula: len(iso) for formula,iso in self.items()}

    def to_frame(self, ions, return_arrays=False):
        """Represent all of isotopic patterns as pd.DataFrame. """
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
