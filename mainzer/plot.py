import numpy as np


def plot_spectrum(mz, intensity, idx=None, show=True, **kwds):
    import matplotlib.pyplot as plt
    plt.stem(mz, intensity, markerfmt=' ', use_line_collection=True, **kwds)
    if idx is not None:
        plt.scatter(mz[idx], intensity[idx], c='black')
    if show:
        plt.show()


def plot_mz_diffs(mz, show=True):
    import matplotlib.pyplot as plt
    plt.plot(mz[:-1], np.diff(mz))
    if show:
        plt.show()

def segments(X, Y):
    ox = []
    oy = []
    for x, y in zip(X, Y):
        ox.append(x); oy.append(0.0)
        ox.append(x); oy.append(y)
        ox.append(x); oy.append(0.0)
    return np.array(ox), np.array(oy)


def get_layout(scheme='inversed'):
    import plotly.graph_objs as go
    if scheme== 'inversed':
        font_c = 'white'
        back_c = 'black'
    else:
        font_c = 'black'
        back_c = 'white'
    layout = go.Layout(font = dict(color = font_c),
                       yaxis = dict(title = "Intensity",
                                    color = font_c,
                                    showline = False,
                                    zeroline = False),
                       xaxis = dict(title = "mass/charge [Th]",
                                    color = font_c),
                       plot_bgcolor = back_c,
                       paper_bgcolor= back_c)
    return layout


def plotly_spectrum(spectrum,
                    path="spectrum.html",
                    sticks=False,
                    color_scheme="inversed",
                    webgl=True,
                    show=True):
    """Make a plotly plot of the fittings.

    The plot overlays scatterplot of fitted intensities
    over the bars corresponding to total intensities in peak groups.

    Parameters
    ==========
    mz : numpy.array
        Mass to charge ratios recorded in the spectrum.
    intensity : numpy.array
        Intensities recorded in the spectrum.
    path : str
        Where to save the plot.
        The file should have 'html' extension to avoid warnings.
    webgl : boolean
        Should we use WebGL? It's quicker to use: try it!
    show : boolean
        Open the browser with the output html file?
    """
    import plotly
    import plotly.graph_objs as go
    scatter = go.Scattergl if webgl else go.Scatter
    mz, intensity = zip(*spectrum.confs)
    lines  = scatter(x=mz, y=intensity,
                     line=dict(color="orange"),
                     name='Observed')
    data = [lines]
    if sticks:
        sx, sy = segments(mz, intensity)
        points = scatter(x=sx, y=sy,
                         line = dict(color='blue'),
                         name='Hidden')
        data.append(points)
    layout = get_layout(color_scheme)
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=path, auto_open=show)


def spec_envelope(spectrum, name=None, line={'color':'orange'}, webgl=True, **kwds):
    import plotly.graph_objs as go
    scatter = go.Scattergl if webgl else go.Scatter
    mz, intensity = zip(*spectrum.confs)
    return scatter(x=mz, y=intensity, line=line, name=name, **kwds)

def spec_sticks(spectrum, name=None, line={'color':'white'}, webgl=True, **kwds):
    import plotly.graph_objs as go
    scatter = go.Scattergl if webgl else go.Scatter
    mz, intensity = zip(*spectrum.confs)
    sx, sy = segments(mz, intensity)
    return scatter(x=sx, y=sy, line=line, name=name, **kwds)

def plotly_deconv(data,
                  path="deconv.html",
                  color_scheme="inversed",
                  webgl=True,
                  show=True):
    import plotly
    import plotly.graph_objs as go
    layout = get_layout(color_scheme)
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=path, auto_open=show)

def plot_histogram1D(data, bins=100, show=True):
    plt.hist(data, bins=bins)
    if show:
        plt.show()

def plot_hexbin(x, y, show=True, **kwds):
    plt.hexbin(x,y,**kwds)
    if show:
        plt.show()

def plot_fit(centroids_df, show=True):
    import matplotlib.pyplot as plt
    plot_spectrum(centroids_df.mz_apex, centroids_df.integrated_intensity, show=False)
    centroids_df_short = centroids_df.query("chimeric_intensity_in_centroid > 0")
    plt.scatter(centroids_df_short.mz_apex, centroids_df_short.chimeric_intensity_in_centroid, c="red")
    if show:
        plt.show()

