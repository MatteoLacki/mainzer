import pandas as pd


def iter_charge_sequence_groups(x, min_len=1):
    x = iter(x)
    prev = next(x)
    cluster = [0]
    for i, curr in enumerate(x, 1):
        if curr == prev + 1:
            cluster.append(i)
        else:
            if len(cluster) >= min_len:
                yield cluster
            cluster = [i]
        prev = curr
    if len(cluster) >= min_len:
        yield cluster


def charge_sequence_filter(ions, min_charge_sequence_length=3):
    ions = ions.copy()
    ions.sort_values(["name","charge"], inplace=True)
    assert all(col in ions.columns for col in ("name","charge")), \
        "Columns should contain 'name' and 'charge'."
    
    sequence_charge_groups = []
    i = 0
    for name, ion in ions.groupby("name"):
        for idx in iter_charge_sequence_groups(ion.charge, min_charge_sequence_length):
            sequence_charge_group = ion.iloc[idx].copy()
            sequence_charge_group["group_id"] = i
            sequence_charge_groups.append(sequence_charge_group)
            i += 1
    return pd.concat(sequence_charge_groups)
