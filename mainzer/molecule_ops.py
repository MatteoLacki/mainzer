import collections
import itertools
import pandas as pd

from .formulas import formula2counter, counter2formula



def iter_mers(molecules, degree=2):
    """Turn molecules into mers of molecules of a given degree.

    Args:
        molecules (dict): maps name to chemical formula (string to string).
        degree (int): How many molecules can merge?

    Yields:
        tuple: pair consisting of a name and a formula that adds the formulas of the component molecules.
    """
    if degree == 0:
        yield from molecules.items()
    else:
        for comb in itertools.combinations_with_replacement(molecules.items(), degree):
            names, formulas = zip(*comb)
             
            cnt = collections.Counter()
            for f in formulas:
                cnt += formula2counter(f)
            name = "_".join(names)
            yield name, counter2formula(cnt)


def iter_mered_molecules(molecules, max_degree):
    """Yield all mered molecule up to a certain degree of merging.

    Args:
        molecules (dict): maps name to chemical formula (string to string).
        max_degree (int): Maximal number of molecules to merge.

    Yields:
        tuple: pair consisting of a name and a formula that adds the formulas of the component molecules.
    """
    for degree in range(max_degree+1):
        yield from iter_mers(molecules, degree)


def crosslink(*args, separator='__'):
    """Make crosslinks between formulas. 
    
    Arguments:
        *args: dictionaries of molecules to crosslink.

    Yields:
        tuple: compound name of crosslink molecule and its formula.
    """
    for prod in itertools.product(*[d.items() for d in args]):
        names, formulas = zip(*prod)
        name = separator.join(names)
        final_formula = collections.Counter()
        for formula in formulas:
            final_formula += formula2counter(formula)
        final_formula = counter2formula(final_formula)
        yield name, final_formula


def molecules2df(molecules, charges=None):
    res = pd.DataFrame(molecules.items(), columns=("name", "formula"))
    return res if charges is None else res.merge(pd.Series(charges, name="charge"), how="cross")
