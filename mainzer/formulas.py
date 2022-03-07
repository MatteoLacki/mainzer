import collections
import re

import aa2atom


def formula2counter(formula: str) -> collections.Counter:
    parsed = re.findall('([A-Z][a-z]*)(-?[0-9]*)', formula)
    cnt = collections.Counter()
    for e, n in parsed:
        n = int(n) if n else 1
        cnt[e] += n
    return cnt


def counter2formula(cnt: collections.Counter, sep: str="") -> str:
    return "".join(f"{el}{sep}{c}" for el,c in cnt.items())


def aa2formula(amino_acid_sequence: str) -> str:
    return counter2formula(aa2atom.aa2atom(amino_acid_sequence))


def add_formulas(formulas_0: str, formulas_1: str):
    formulas_0_cnt = formula2counter(formulas_0)
    formulas_1_cnt = formula2counter(formulas_1)
    return counter2formula(formulas_0_cnt + formulas_1_cnt)