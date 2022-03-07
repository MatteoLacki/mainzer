"""This is like a wrapper to initialize the formulas."""
import IsoSpecPy
import numpy as np
import pandas as pd
from typing import Dict, Tuple, Iterable, List, Union

from dataclasses import dataclass



@dataclass
class IsotopicCalculator:
    coverage: float
    bin_size: float
    get_minimal_pset: bool=True
    PROTON_MASS: float=1.007276
    
    def __call__(self, formula: str) -> IsoSpecPy.IsoDistribution:
        envelope = IsoSpecPy.IsoTotalProb(
            prob_to_cover=self.coverage,
            formula=formula,
            get_minimal_pset=self.get_minimal_pset
        )
        if self.bin_size is not None:
            envelope = envelope.binned(self.bin_size)
        return envelope

    def masses_probs(self, formula: str) -> Tuple[np.array]:
        envelope = self(formula)
        masses = envelope.np_masses()
        probs = envelope.np_probs()
        return masses, probs

    def spectrum(self, formula: str, charge: int=1) -> Tuple[np.array]:
        envelope = self(formula)
        masses = envelope.np_masses()
        probs = envelope.np_probs()
        return (masses / charge + self.PROTON_MASS, probs)

    def masses(self, formula: str) -> np.array:
        return self(formula).np_masses()

    def mzs(self, formula: str, charge: int=1) -> np.array:
        return self.masses(formula) / charge + self.PROTON_MASS

    def probs(self, formula: str) -> np.array:
        return self(formula).np_probs()    

    def iter_envelope_summaries(
        self, 
        formulas: Iterable[str]
    ) -> Iterable[Dict]:
        for formula in formulas:
            envelope = self(formula)
            masses = envelope.np_masses()
            probs = envelope.np_probs()
            top_prob_idx = np.argmax(envelope.np_probs())
            yield {"formula": formula,
                   "envelope_size": len(envelope),
                   "envelope_total_prob": probs.sum(),
                   "min_mass": masses.min(),
                   "max_mass": masses.max(),
                   "top_prob_mass": masses[top_prob_idx]}

    def summary_df(
        self,
        formulas: Iterable[str]
    ) -> pd.DataFrame:
        return pd.DataFrame(self.iter_envelope_summaries(formulas))

    def ions_summary(
        self,
        ions_df: pd.DataFrame,
        copy_ions_df: bool=True
    ) -> pd.DataFrame:
        assert all(col in ions_df.columns for col in ("formula", "charge"))
        if copy_ions_df:
            ions_df = ions_df.copy()
        res = self.summary_df(ions_df.formula.unique())
        ions_df = ions_df.merge(res)
        ions_df["envelope_min_mz"] = ions_df.min_mass / ions_df.charge + self.PROTON_MASS
        ions_df["envelope_max_mz"] = ions_df.max_mass / ions_df.charge + self.PROTON_MASS
        ions_df["envelope_top_prob_mz"] = ions_df.top_prob_mass / ions_df.charge + self.PROTON_MASS
        ions_df.drop(["min_mass","max_mass","top_prob_mass"], axis=1, inplace=True)
        return ions_df

    # def sizes(self) -> :
    #     return {formula: len(iso) for formula,iso in self.items()}

    def to_frame(
        self,
        formula: str,
        charge: int
    ) -> pd.DataFrame:
        mz, prob = self.spectrum(formula, charge)
        return pd.DataFrame(
            {
                "formula":formula,
                "charge":charge,
                "mz":mz,
                "prob":prob
            }
        )

    def big_frame(
        self,
        ions: List[Tuple[str,int]],
        return_arrays: bool=False
    ) -> Union[Tuple[np.array], pd.DataFrame]:
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
            res = pd.DataFrame(
                {
                    "formula": np.repeat(
                            ions.get_level_values("formula"),
                            theory_peaks_cnts
                        ),
                    "charge": np.repeat(ions.get_level_values("charge"), theory_peaks_cnts),
                    "isospec_mz":mzs,
                    "isospec_prob":probs})
            return res
