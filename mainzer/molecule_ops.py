import collections
import itertools
import pandas as pd

from .formulas import formula2counter, counter2formula, aa2formula, add_formulas



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
            print(counter2formula(collections.Counter(names), "_"))

            cnt = collections.Counter()
            for f in formulas:
                cnt += formula2counter(f)
            # name = "_".join(names)
            name = f"{names[0]}_{len(names)}"
            yield name, counter2formula(cnt)


def iter_mered_molecules(molecules, min_degree, max_degree):
    for degree in range(min_degree, max_degree+1):
        yield from iter_mers(molecules, degree)


def iter_mers2(molecules, degree=2, cnt_sep=" ", mer_sep=" + "):
    for comb in itertools.combinations_with_replacement(molecules.items(), degree):
        formula_cnt = collections.Counter()
        names_cnt = collections.Counter()
        for mol_name, mol_formula in comb:
            formula_cnt += formula2counter(mol_formula)
            names_cnt[mol_name] += 1
        name = mer_sep.join(f"{cnt}{cnt_sep}{name}" if cnt>1 else name for name, cnt in names_cnt.items())
        formula = counter2formula(formula_cnt)
        yield name, formula

def iter_mered_molecules2(molecules, degrees=[1]):
    for degree in degrees:
        yield from iter_mers2(molecules, degree)




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


def mered_ions(molecules, mer_degrees, charges):
    return molecules2df(dict(iter_mered_molecules2(molecules, mer_degrees)), charges)

def mered_proteins(base_proteins, only_heteromers=False):
    homomers = []
    for base_protein in base_proteins.itertuples():
        homomers.append(dict(iter_mered_molecules2({base_protein.name : base_protein.formula }, base_protein.oligomeric_states)))
    heteromers = dict(crosslink(*homomers, separator=" "))
    if only_heteromers:
        return heteromers
    else:
        all_mers = heteromers 
        for h in homomers:
            all_mers.update(h)
        return all_mers


def mered_lipids(base_lipids, min_lipid_mers, max_lipid_mers):
    base_lipids_dicts = dict(zip(base_lipids.name, base_lipids.formula))
    mer_degs = range(min_lipid_mers, max_lipid_mers+1)
    return dict(iter_mered_molecules2(base_lipids_dicts, mer_degs))


def merge_free_lipids_and_promissing_proteins(free_lipids,
                                              promissing_proteins):
    assert all(col in free_lipids.columns for col in ("name","formula"))
    assert all(col in promissing_proteins.columns for col in ("name","formula","charge"))
    res = promissing_proteins[["name","formula","charge"]]
    prot_col_names = [f"prot_{col}" for col in res.columns]    
    res.columns = prot_col_names
    res = pd.merge(res.assign(key=0),
                             free_lipids.assign(key=0),
                             how="outer", on="key").drop("key", axis=1)
    
    res.name = res.prot_name + " + " +  res.name
    res.formula = [add_formulas(f_prot, f_lip) for f_prot, f_lip in zip(res.prot_formula, res.formula)]
    res = res.drop(['prot_formula', 'prot_name'], axis=1)
    res.rename(columns={"prot_charge":"charge"}, inplace=True)
    res = res[["name","formula","charge"]]
    return res