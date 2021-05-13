import aa2atom
import pandas as pd

from .molecule_ops import iter_mered_molecules, crosslink, molecules2df


fasta_2_formula = lambda aa: aa2atom.atom2str(aa2atom.aa2atom(aa))


def get_lipido_ions(molecules,
                    max_protein_mers,
                    max_lipid_mers,
                    min_lipid_charge,
                    max_lipid_charge,
                    min_protein_charge,
                    max_protein_charge,
                    **kwargs):
    """Get a set of ions for lipido.

    Arguments:
        molecules (pd.DataFrame): A data frame with columns "group", "name", and "sequence_or_formula".
        max_protein_mers (int):
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
    protein_dicts = [aa2atom.aa2atom(fasta) for fasta in proteins.sequence_or_formula]

    pmers_lst = []
    # Only allowing straight mers, not complexes of different proteins
    for name, protein_dict in zip(proteins.name, protein_dicts):
        formula = protein_dict
        for n in range(1, max_protein_mers+1):
            rname = name if n == 1 else str(n) + "meric-" + name
            pmers_lst.append(molecules2df({rname: aa2atom.atom2str(formula)}, range(n*min_protein_charge, n*max_protein_charge+1)))
            formula = formula + protein_dict

    protein_formulas = [fasta_2_formula(fasta) for fasta in proteins.sequence_or_formula]
    proteins = dict(zip(proteins.name, protein_formulas))
    lipids = molecules.query("group == 'lipid'")
    lipids = dict(zip(lipids.name, lipids.sequence_or_formula))
    lipid_mers = dict(iter_mered_molecules(lipids, max_lipid_mers))
    lipid_protein_mers = dict(crosslink(lipid_mers, proteins))


    ions = pd.concat([molecules2df(lipid_mers, range(min_lipid_charge, max_lipid_charge+1)),
                      molecules2df(lipid_protein_mers, range(min_protein_charge, max_protein_charge+1))] +\
                      pmers_lst,
                      ignore_index=True)
    return ions



def get_roepstorrf_scheme_fragments(proteins, fragments="by"):
    """Generate fragments.

    Arguments:
        proteins (pd.DataFrame): A data frame with columns "name" and "sequence".
        fragments (str): Types of fragments to generate.
    """

    for fragment_type in fragments:
        pass
    pass
