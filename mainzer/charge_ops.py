import numpy as np # needed only for array and arg_sort


def cluster_sorted_charges(sorted_charges_array, min_charge_sequence=3):
    charge_prev = sorted_charges_array[0]
    assignments = np.zeros(len(sorted_charges_array), dtype=int)
    cluster_No  = 1
    cluster = []
    for i, charge in enumerate(sorted_charges_array):
        if charge != charge_prev + 1:
            if len(cluster) >= min_charge_sequence:
                assignments[cluster] = cluster_No
                cluster_No += 1
            cluster = []
        cluster.append(i)
        charge_prev = charge
    if len(cluster) >= min_charge_sequence:
        assignments[cluster] = cluster_No
    return assignments


def cluster_charges(charges, min_charge_sequence=3):
    charges = np.array(charges)
    order = np.argsort(charges)
    sorted_charges_array = charges[order]
    assignments = cluster_sorted_charges(
        sorted_charges_array, 
        min_charge_sequence
    )
    return assignments[np.argsort(order)]


def test_charge_clustering():
    charges = [6, 1,2,3,4,  -1,  8,9,10, 12, 16,15, 13, 18, 17, 20]
    assert np.all(cluster_charges(charges) == np.array([0, 1,1,1,1, 0, 2,2,2, 0, 3,3, 0, 3,3, 0]))


def test_charge_clustering_ions():
    import pandas as pd # only needed for the test
    
    formulas_charges = pd.DataFrame({
        "formula": ['a','a','a','b', 'a','b','b','b','b','b','b', 'c','c'],
        "charge" : [ 1,  2,  3,  7,   6,  2,  3,  6,  1,  8,  10,  1,  2]}
    )
    formulas_charges["assignment"] = formulas_charges.groupby("formula").charge.transform(cluster_charges)

    assert np.all(formulas_charges.assignment.values == np.array([1,1,1,2,0,1,1,2,1,2,0,0,0]))