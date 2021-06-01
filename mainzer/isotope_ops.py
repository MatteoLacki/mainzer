"""This is like a wrapper to initialize the formulas."""
import numpy as np
import pandas as pd

import IsoSpecPy as isospec


class IsotopicCalculator(object):
    def __init__(self,
                 coverage=.95,
                 bin_size=.1,
                 get_minimal_pset=True,
                 PROTON_MASS=1.007276):
        self.coverage = coverage
        self.bin_size = bin_size
        self.get_minimal_pset = get_minimal_pset
        self.PROTON_MASS = PROTON_MASS
        

    def __call__(self, formula):
        envelope = isospec.IsoTotalProb(prob_to_cover=self.coverage,
                                        formula=formula,
                                        get_minimal_pset=self.get_minimal_pset)
        if self.bin_size is not None:
            envelope = envelope.binned(self.bin_size)
        return envelope

    def spectrum(self, formula, charge=1):
        envelope = self(formula)
        masses = envelope.np_masses()
        probs = envelope.np_probs()
        if charge > 1:
            return (masses / charge + self.PROTON_MASS, probs)
        else:
            return masses, probs

    def masses(self, formula):
        return self(formula).np_masses()

    def mzs(self, formula, charge=1):
        return self.masses(formula) / charge + self.PROTON_MASS

    def probs(self, formula):
        return self(formula).np_probs()    

    def iter_envelope_summaries(self, formulas):
        for formula in formulas:
            envelope = self(formula)
            masses = envelope.np_masses()
            probs = envelope.np_probs()
            yield {"formula": formula,
                   "envelope_size": len(envelope),
                   "envelope_total_prob": probs.sum(),
                   "min_mass": masses.min(),
                   "max_mass": masses.max()}

    def summary_df(self, formulas):
        return pd.DataFrame(self.iter_envelope_summaries(formulas))

    def ions_summary(self, ions_df, copy_ions_df=True):
        assert all(col in ions_df.columns for col in ("formula", "charge"))
        if copy_ions_df:
            ions_df = ions_df.copy()
        res = self.summary_df(ions_df.formula.unique())
        ions_df = ions_df.merge(res)
        ions_df["envelope_min_mz"] = ions_df.min_mass / ions_df.charge + self.PROTON_MASS
        ions_df["envelope_max_mz"] = ions_df.max_mass / ions_df.charge + self.PROTON_MASS
        ions_df.drop(['min_mass','max_mass'], axis=1, inplace=True)
        return ions_df

    def sizes(self):
        return {formula: len(iso) for formula,iso in self.items()}

    def to_frame(self, formula, charge):        
        mz, prob = self.spectrum(formula, charge)
        return pd.DataFrame({"formula":formula, "charge":charge, "mz":mz, "prob":prob})

    