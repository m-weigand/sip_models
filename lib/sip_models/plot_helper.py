# -*- coding: utf-8 -*-
""" Import all necessary matplotlib modules and set default options To use this
module, import it as:

import sip_models.plot_helper
plt, mpl = sip_models.plot_helper.setup()
"""
import os
import sys
import platform


def setup():
    """import the matplotlib modules and set the style

    Returns
    -------
    plt: pylab
        imported pylab module
    mpl: matplotlib module
        imported matplotlib module

    """
    # Latex support can be activated using an environment variable, otherwise
    # the default settings are:
    # - for windows: off
    # - else: on
    use_latex = False
    if('DD_USE_LATEX' in os.environ):
        if os.environ['DD_USE_LATEX'] == '1':
            use_latex = True
    else:
        if platform.system() == "Windows":
            use_latex = False
        else:
            use_latex = True

    already_loaded = 'matplotlib' in sys.modules

    # just make sure we can access matplotlib as mpl
    import matplotlib as mpl

    if not already_loaded:
        mpl.use('Agg')

    import matplotlib.pyplot as plt

    plt.style.use('seaborn')

    # general settings
    mpl.rcParams['font.size'] = 7.0
    mpl.rcParams['axes.labelsize'] = 7.0
    mpl.rcParams['xtick.labelsize'] = 7.0
    mpl.rcParams['ytick.labelsize'] = 7.0
    mpl.rcParams["lines.linewidth"] = 1.5
    mpl.rcParams["lines.markeredgewidth"] = 3.0
    mpl.rcParams["lines.markersize"] = 3.0
    # mpl.rcParams['font.sans-serif'] = 'Droid Sans'

    # mpl.rcParams['font.family'] = 'Open Sans'
    # mpl.rcParams['font.weight'] = 400
    mpl.rcParams['mathtext.default'] = 'regular'

    # mpl.rcParams['font.family'] = 'Droid Sans'

    if use_latex:
        mpl.rcParams['text.usetex'] = True

        mpl.rc(
            'text.latex',
            preamble=''.join((
                # r'\usepackage{droidsans}',
                # r'\usepackage[T1]{fontenc} ',
                r'\usepackage{sfmath} \renewcommand{\rmfamily}{\sffamily}',
                r'\renewcommand\familydefault{\sfdefault} ',
                # r'\usepackage{mathastext} '
            ))
        )
    else:
        mpl.rcParams['text.usetex'] = False
    import mpl_toolkits.axes_grid1 as axes_grid1
    axes_grid1
    return plt, mpl


def mpl_get_cb_bound_next_to_plot(ax):
    """
    Return the coordinates for a colorbar axes next to the provided axes
    object. Take into account the changes of the axes due to aspect ratio
    settings.

    Parts of this code are taken from the transforms.py file from matplotlib

    Important: Use only AFTER fig.subplots_adjust(...)

    Use as:
    =======
    cb_pos = mpl_get_cb_bound_next_to_plot(fig.axes[15])
    ax1 = fig.add_axes(cb_pos, frame_on=True)

    cmap = mpl.cm.jet_r
    norm = mpl.colors.Normalize(vmin=float(23), vmax=float(33))
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                                    orientation='vertical')
    cb1.locator = mpl.ticker.FixedLocator([23,28,33])
    cb1.update_ticks()
    cb1.ax.artists.remove(cb1.outline)    # remove framei
    """
    position = ax.get_position()

    figW, figH = ax.get_figure().get_size_inches()
    fig_aspect = figH / figW
    box_aspect = ax.get_data_ratio()
    pb = position.frozen()
    pb1 = pb.shrunk_to_aspect(box_aspect, pb, fig_aspect).bounds

    ax_size = ax.get_position().bounds

    xdiff = (ax_size[2] - pb1[2]) / 2
    ydiff = (ax_size[3] - pb1[3]) / 2

    # the colorbar is set to 0.01 width
    sizes = [ax_size[0] + xdiff + ax_size[2] + 0.01,
             ax_size[1] + ydiff, 0.01, pb1[3]]

    return sizes


def mpl_get_cb_bound_below_plot(ax):
    """
    Return the coordinates for a colorbar axes below the provided axes object.
    Take into account the changes of the axes due to aspect ratio settings.

    Parts of this code are taken from the transforms.py file from matplotlib

    Important: Use only AFTER fig.subplots_adjust(...)

    Use as:
    =======
    """
    position = ax.get_position()

    figW, figH = ax.get_figure().get_size_inches()
    fig_aspect = figH / figW
    box_aspect = ax.get_data_ratio()
    pb = position.frozen()
    pb1 = pb.shrunk_to_aspect(box_aspect, pb, fig_aspect).bounds

    ax_size = ax.get_position().bounds

    # the colorbar is set to 0.01 width
    sizes = [ax_size[0], ax_size[1] - 0.14, pb1[2], 0.03]

    return sizes
