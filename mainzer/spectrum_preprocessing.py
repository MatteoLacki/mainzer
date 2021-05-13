from .models import fit_polynomial


def get_mz2dmz(mzs, resolutions):
    mz2res, coefs = fit_polynomial(mzs, resolutions, 1)
    def mz2dmz(mz):
        return mz/mz2res(mz)
    return mz2dmz, mz2res, coefs


def mz2dmz_ppm(mz, ppm=10):
    return mz*ppm/10**6


def trim_low_intense_peaks(mz, intensity, cut_off):
    """Trim peaks that are not intense enough.
    
    Arguments:
        mz (np.array):
        intensity (np.array):
        cut_off (float):
    """
    return mz[intensity >= cut_off], intensity[intensity >= cut_off]
