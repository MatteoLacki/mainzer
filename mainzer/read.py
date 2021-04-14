import pathlib
import sys
import pandas as pd



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


def read(path):
    path = pathlib.Path(path)
    extension = path.suffix.lower()
    if extension in ('.mzml','.mzxml'):
        return mzml(path)
    elif extension in ('.csv'):    
        return csv(path)
    else:
        raise NotImplementedError(f"We don't know how to parse {path}")
