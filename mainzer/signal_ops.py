import numpy as np
import pandas as pd


def iter_rel_minima(intensity):
    I = iter(intensity)
    i__ = next(I)
    _i_ = next(I)
    if i__ < _i_:
        yield 0
    for idx, __i in enumerate(I,2):
        if i__ > _i_ <= __i or i__ >= _i_ < __i :
            yield idx-1
        _i_, i__ = __i, _i_
    if __i < _i_:
        yield idx 


def iter_peaks(intensity):
    """Divide the peaks based on bitonicity. """
    I = iter(intensity)
    i__ = next(I)
    _i_ = next(I)
    peakNo = 0
    yield peakNo
    if i__ < _i_:
        peakNo += 1
    for idx, __i in enumerate(I,2):
        yield peakNo
        if i__ > _i_ <= __i or i__ >= _i_ < __i :
            peakNo += 1
        _i_, i__ = __i, _i_
    if __i < _i_:
        peakNo += 1
    yield peakNo


def rel_minima(intensity):
    return np.array(list(iter_rel_minima(intensity)))

#TODO: there should be a distance parameter?
def peak_clusters(intensity):
    return np.array(list(iter_peaks(intensity)))


def iter_borders(intensity):
    it = iter_rel_minima(intensity)
    i0 = next(it)
    if i0 != 0:
        yield 0
    yield i0
    for i in it:
        yield i
    if i < len(intensity):
        yield len(intensity)


def iter_paired_borders(intensity):
    it = iter_borders(intensity)
    i_ = next(it)
    for _i in it:
        yield i_, _i
        i_ = _i


def iter_grouped_argmax(clusters, Y):
    """Iterate grouped arg-maxes.

    Args:
        clusters (iterable): a non-decreasing sequence of integers
    Yields:
        int: position of maxima.
    """
    max_y = 0
    max_y_idx = 0
    prev_cl = -1
    for i, (cl, y) in enumerate(zip(clusters, Y)):
        if cl == prev_cl + 2: 
            yield max_y_idx
            max_y = y
            max_y_idx = i
            prev_cl = cl-1
        else:
            if y > max_y:
                max_y = max(y, max_y)
                max_y_idx = i
    yield max_y_idx


def conditional_argmax(Y, clusters):
    """Return grouped arg-maxes.

    Args:
        Y (iterable): values for which to calculate conditional arg-maxima.
        clusters (iterable): a non-decreasing sequence of integers
    Returns:
        np.array: indices of maxima in clusters.
    """
    return np.fromiter(iter_grouped_argmax(clusters, Y), int)


def centroid(x, I):
    """Centroid x according to I.

    Args:
        x (list): positions to center.
        I (list): intensities for centroiding.
    Returns:
        tuple: positions of peak apices, their intensity, the intensity of an entire peak, clustering of intensity indices, and indices of apexes. 
    """
    x = np.array(x)
    I = np.array(I)
    assert np.all(np.diff(x) >= 0), "x must be non-decreasing"
    assert np.all(I >= 0), "I must be non-negative."
    clusters = peak_clusters(I)
    max_I_idx = conditional_argmax(I, clusters)
    max_I = I[max_I_idx]
    sum_I = np.bincount(clusters, I)
    x_apex = x[max_I_idx]    
    return x_apex, max_I, sum_I, clusters, max_I_idx


def get_cluster_ends(clusters):
    """Get the indices of beginnings and ends of clusters.

    Args:
        clusters (np.array): nondecreasing sequence of integers that groups consecutive indices into clusters. 

    Returns:
        tuple: start and end indices for clusters.
    """
    prev_cl = clusters[0]
    left_idxs = [0]
    right_idxs = []
    for i, CL in enumerate(clusters):
        if CL != prev_cl:
            right_idxs.append(i-1)
            left_idxs.append(i)
        prev_cl = CL
    right_idxs.append(i)
    start_idxs = np.array(left_idxs)
    end_idxs = np.array(right_idxs)
    return start_idxs, end_idxs


def weighted_means_and_std_devs(mz, intensity, mz_clusters, I_sum):
    wMZs = np.bincount(mz_clusters, mz * intensity)[I_sum > 0] / I_sum[I_sum > 0]
    wMS2z = np.bincount(mz_clusters, mz**2 * intensity)[I_sum > 0] / I_sum[I_sum > 0]
    stds = np.sqrt(wMS2z - wMZs**2)
    return wMZs, stds


def estimate_mz2std(mz, intensity):
    mz_apex, I_max, I_sum, mz_clusters, I_max_idx = centroid(mz, intensity)
    wMZs, stds = weighted_means_and_std_devs(mz, intensity, mz_clusters, I_sum)
    mz2std = np.poly1d(np.polyfit(wMZs, stds, deg=3))
    return mz2std


def cluster_spectrum(mz, intensity):
    mz_apex, I_max, I_sum, mz_clusters, I_max_idx = centroid(mz, intensity)
    clustered_spectrum = pd.DataFrame({"mz_apex":mz_apex, "I_max":I_max, "I_sum": I_sum})
    clustered_spectrum["left_idx"], clustered_spectrum["right_idx"] = get_cluster_ends(mz_clusters)
    clustered_spectrum["left_mz"] = mz[clustered_spectrum.left_idx]
    clustered_spectrum["right_mz"] = mz[clustered_spectrum.right_idx]
    clustered_spectrum.index.name = 'cluster'
    return clustered_spectrum
