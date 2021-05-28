import pathlib
import sys
import pandas as pd

from .formulas import aa2formula



def mzml(path, verbose=False):
    import pyteomics.mzml
    import pyteomics.mzxml
    path = pathlib.Path(path)
    ext = path.suffix.lower()
    if ext in ('.mzml','.mzxml'):
        if ext == '.mzml':
            read = pyteomics.mzml.read
        else:
            read = pyteomics.mzxml.read
        data = next(read(str(path)))
        if verbose:
            print(data)
        return data['m/z array'], data['intensity array']
    else:
        raise NotImplementedError(f"We don't know how to parse {path}")


def csv(path):
    spectrum = pd.read_csv(path)
    return spectrum.iloc[:,0], spectrum.iloc[:,1]


def read_spectrum(path):
    path = pathlib.Path(path)
    extension = path.suffix.lower()
    if extension in ('.mzml','.mzxml'):
        return mzml(path)
    elif extension in ('.csv'):    
        return csv(path)
    else:
        raise NotImplementedError(f"We don't know how to parse {path}")


def read_molecules_for_lipido(path):
    molecules = pd.read_csv(path)
    molecules["formula"] = molecules["sequence_or_formula"]
    molecules.loc[molecules.group == "protein", "formula"] = molecules[molecules.group == "protein"].formula.apply(aa2formula)
    proteins = molecules[molecules.group == "protein"]
    proteins = dict(zip(proteins.name, proteins.formula))
    lipids = molecules[molecules.group == "lipid"]
    lipids = dict(zip(lipids.name, lipids.formula))
    return proteins, lipids


def parse_oligomeric_state(state, sequence_sep="-"):
    res = set([])
    try:
        for term in state.replace(" ","").split(","):
            if sequence_sep in term:
                start, end = term.split(sequence_sep)
                res.update(range(int(start), int(end)+1))
            else:
                res.add(int(term))
    except ValueError as ve:
        print(repr(ve))
        print("Allowed are expressions with integers and 'integer:another_integer'.")
        raise ValueError
    res = list(res)
    res.sort()
    return res

# path = "/home/matteo/Projects/och_Kallol/mainzer/test_data/base_proteins.csv"
def read_base_proteins(path):
    base_proteins = pd.read_csv(path)
    base_proteins["formula"] = base_proteins["amino_acidic_sequence"].apply(aa2formula)
    base_proteins["oligomeric_states"] = base_proteins.oligomeric_states.apply(parse_oligomeric_state)
    return base_proteins

# path = "/home/matteo/Projects/och_Kallol/mainzer/test_data/base_lipids.csv"
def read_base_lipids(path):
    base_lipids = pd.read_csv(path)
    return base_lipids