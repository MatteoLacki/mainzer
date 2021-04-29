import matplotlib.pyplot as plt

def plot_spectrum(mz, intensity, idx=None, show=True, **kwds):
    import matplotlib.pyplot as plt
    plt.stem(mz, intensity, markerfmt=' ', use_line_collection=True, **kwds)
    if idx is not None:
        plt.scatter(mz[idx], intensity[idx], c='black')
    if show:
        plt.show()

plot_spectrum(mz, intensity)