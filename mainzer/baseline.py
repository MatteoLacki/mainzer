import numpy as np
from typing import Tuple


def strip_baseline(
    mzs: np.array,
    intensities: np.array,
    baseline: float
) -> Tuple[np.array]:
    if len(mzs) == 0:
        return [], []

    idx = 0
    forward_idx = 0
    lookahead = 10.0
    dens = baseline / lookahead

    res_mzs = []
    res_int = []


    while idx < len(mzs):
        while forward_idx < len(mzs) and mzs[idx] + lookahead > mzs[forward_idx]:
            forward_idx += 1

        if forward_idx < len(mzs):
            try:
                dens = baseline * lookahead / (forward_idx - idx)
            except ZeroDivisionError:
                pass


        if intensities[idx] >= dens:
            res_mzs.append(mzs[idx])
            res_int.append(intensities[idx])

        idx += 1
        forward_idx = max(idx, forward_idx)

    return np.array(res_mzs), np.array(res_int)

