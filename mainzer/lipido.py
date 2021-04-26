import pandas as pd
from aa2atom import aa2atom, atom2str
aa2str = lambda aa: atom2str(aa2atom(aa))

from .molecule_ops import iter_mered_molecules, crosslink, molecules2df

def get_lipido_ions(molecules,
                    max_lipid_mers,
                    min_lipid_charge,
                    max_lipid_charge,
                    min_protein_charge,
                    max_protein_charge,
                    **kwargs):
    """Get a set of ions for lipido.

    Arguments:
        molecules (pd.DataFrame): A data frame with columns "group", "name", and "sequence_or_formula".
        max_lipid_mers (int): 
        min_lipid_charge (int):
        max_lipid_charge (int):
        min_protein_charge (int):
        max_protein_charge (int):
        **kwargs: other key-world arguments.
    Returns:
        pd.DataFrame: A data.frame with columns "name", "formula", and "charge": the general input for intensity estimation algorithms.
    """
    assert all(col in molecules.columns for col in ("group","name","sequence_or_formula")), "The csv with molecules should have columns 'group', 'name', and 'sequence'."
    proteins = molecules.query("group == 'protein'")
    protein_formulas = []
    for seq in proteins.sequence_or_formula:
        formula = aa2str(seq)
        protein_formulas.append(formula)
        # except Exception:
        #     protein_formulas.append(seq)

    proteins = dict(zip(proteins.name, protein_formulas))
    lipids = molecules.query("group == 'lipid'")
    lipids = dict(zip(lipids.name, lipids.sequence_or_formula))

    lipid_mers = dict(iter_mered_molecules(lipids, max_lipid_mers))
    lipid_protein_mers = dict(crosslink(lipid_mers, proteins))

    ions = pd.concat([molecules2df(lipid_mers, range(min_lipid_charge, max_lipid_charge+1)),
                      molecules2df(proteins, range(min_protein_charge, max_protein_charge+1)),
                      molecules2df(lipid_protein_mers, range(min_protein_charge, max_protein_charge+1))],
                      ignore_index=True)
    return ions
