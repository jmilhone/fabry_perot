
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable


def tableau20_colors():
    """Good set of colors for plotting

    Returns:
        list: list of r,g,b color values
    """
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

    return tableau20


def my_hist(ax, data, bins=None, horizontal=False):
    """Custom histogram function for plotting multinest output

    Args:
        ax (matplotlib.axes.Axes): axes to plot with
        data (np.ndarray): array of values to histogram
        bins (int, optional): number of bins in histogram
        horizontal (bool): histogram is horizontal if True

    Returns:
        np.ndarray: bin edges for the histogram
    """
    if bins is not None:
        hist, bins = np.histogram(data, density=True, bins=bins)
    else:
        hist, bins = np.histogram(data, density=True, bins='auto')

    hist *= 100.0

    bw = bins[1] - bins[0]

    if horizontal:
        ax.barh(bins[0:-1], hist * bw, height=bw)#, color='dimgray')  # , alpha=0.5)
        if data.max() > 1000:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        else:
            ax.get_yaxis().get_major_formatter().set_scientific(True)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
    else:
        ax.bar(bins[0:-1], hist * bw, width=bw)#, color='dimgray')  # , alpha=0.5)
        if data.max() > 1000:
            # I don't think this works
            # ax.get_xaxis().get_major_formatter().set_scientific(True)
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        else:
            ax.get_xaxis().get_major_formatter().set_scientific(True)

        ax.get_xaxis().get_major_formatter().set_useOffset(False)
    return bins


def my_hist2d(ax, data1, data2, bins=None, z=30):
    """Custom 2d histrogram option for plotting correlations

    Args:
        ax (matplotlib.axes.Axes): axes to plot with
        data1 (np.ndarray): array of values to histogram
        data2 (np.ndarray): other array of values
        bins (int, optional): number of bins in histogram
        z (int): number of filled contours

    Returns:
        matplotlib.colorbar.Colorbar: colorbar from the 2d histogram
    """
    if bins is not None:
        hist, xx, yy = np.histogram2d(data1, data2, normed=True, bins=bins)
    else:
        hist, xx, yy = np.histogram2d(data1, data2, normed=True)

    dx = xx[1] - xx[0]
    dy = yy[1] - yy[0]

    im = ax.contourf(yy[0:-1], xx[0:-1], hist * dx * dy, z)
    # ax_divider = make_axes_locatable(ax)
    # cax = ax_divider.append_axes("right", size="7%", pad="2%")
    cb = plt.colorbar(im, ax=ax)

    if data1.max() > 1000:
        # ax.get_yaxis().get_major_formatter().set_scientific(True)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    else:
        ax.get_yaxis().get_major_formatter().set_scientific(False)
    if data2.max() > 1000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        # ax.get_xaxis().get_major_formatter().set_scientific(True)
    else:
        ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    return cb


class ClickBuilder:
    """Little class for interactive plot features

    Attributes:
        x (list): x values from clicking
        y (list): y values from clicking
        fig (matplotlib.figure.Figure): matplotlib figure to get clicks from
        cid (int): id for the matplotlib callback
    """

    def __init__(self, fig):
        self.x = []
        self.y = []
        self.fig = fig
        self.cid = fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        """Callback function for registering matplotlib clicks. X and Y values are stored in self.x and self.y
        
        Args:
            event (matplotlib.backend_bases.LocationEvent): matplotlib location event
        """
        self.x.append(event.xdata)
        self.y.append(event.ydata)


def center_plot(data, x0=None, y0=None):
    """Runs clickable center plot for an image

    Args:
        data (np.ndarray): pixel values for image
        x0 (float, default=None): x center for zoomed image, if
            none will use the center of data array
        y0 (float, default=None): y center for zoomed image, if
            none will use the cneter of data array
    Returns:
        Tuple (float, float): x value of click guess,  y value of click guess
    """
    if x0 is None:
        x0 = data.shape[1] / 2.
    if y0 is None:
        y0 = data.shape[0] / 2.
    dx = int(0.3 * x0)
    dy = int(0.3 * y0)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    cb = axs[0].imshow(data, cmap='gray', origin='lower')
    axs[0].add_patch(Rectangle((x0 - dx, y0 - dy), 2 * dx, 2 * dy, lw=2, linestyle='--', ec='red', fc='none'))
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(cb, cax=cax)
    axs[1].imshow(data, cmap='gray', origin='lower')
    axs[1].set_title('Click for center guess')
    plt.setp(list(axs[1].spines.values()), color='red')
    [i.set_linewidth(5) for i in axs[1].spines.values()]
    [i.set_linestyle('--') for i in axs[1].spines.values()]
    axs[1].set_ylim(y0 - dy, y0 + dy)
    axs[1].set_xlim(x0 - dx, x0 + dy)

    center_guess = ClickBuilder(fig)
    plt.tight_layout()
    plt.show()
    fig.canvas.mpl_disconnect(center_guess.cid)
    return center_guess.x[0], center_guess.y[0]


def ringsum_click(r, sig, title='Click Me!'):
    """Runs clickable ringsum plot

    Args:
        r (np.ndarray): bin values
        sig (np.ndarray): ringsum values
        title (str): instructions for clicker displayed as plot title

    Returns:
        Tuple (list, list): x value(s) of click(s),  y value(s) of click(s)
    """
    fig, axs = plt.subplots(figsize=(10, 6))
    axs.plot(r, sig, lw=2)
    axs.set_title(title, fontsize=22)
    axs.set_xlabel('R (pixels)', fontsize=18)
    axs.set_ylabel('Counts', fontsize=18)
    clicks = ClickBuilder(fig)
    plt.show()
    fig.canvas.mpl_disconnect(clicks.cid)
    return clicks.x, clicks.y


def peak_plot(r, sig, peaks, peaks_sd, orders, fax=None, anspks=None, anspks_sd=None):
    """Plots ringsum data with labeled peaks and orders

    Args:
        r (np.ndarray): r array from ringsum
        sig (np.ndarray): ringsum signal
        peaks (dict): dictionary of peak locations
            the keys are the wavelengths
        peaks_sd (dict): dictionary of peak errors
        orders (dict): dictionary of peak orders
            the keys are the wavelengths
        fax (tuple, default=None): the figure and axis
            handles for the plot, if adding to existing
            plot. Default is None, which will make a
            new set of figure and axes and run plt.show()
        anspks (dict): dictionary with peak answers from multinest
        anspks_sd (dict): dictionary with peak sd answers from multinest
    """
    if fax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig, ax = fax

    colors = tableau20_colors()
    ax.plot(r ** 2, sig, 'o-', color=colors[0], lw=2)
    i = 1
    for key in list(peaks.keys()):
        for j, pk in enumerate(peaks[key]):
            pk_sd = 2. * pk * peaks_sd[key][j] / 2.0  # binwidth / 2.0
            ax.axvspan(pk ** 2 - pk_sd, pk ** 2 + pk_sd, color=colors[i],
                       label='{0}: j={1}'.format(key, orders[key][j]), alpha=0.7)
            i += 1
            if i > len(colors):
                i = 0
            if anspks is not None:
                anspk = anspks[key][j]
                anspk_sd = 2. * anspk * anspks_sd[key][j]
                ax.axvspan(anspk ** 2 - anspk_sd, anspk ** 2 + anspk_sd, color=colors[i],
                           label='{0}: j={1} multinest'.format(key, orders[key][j]), alpha=0.7)
                i += 1
                if i > len(colors):
                    i = 0
    ax.legend(fontsize=14)
    ax.set_xlabel(r'R$^{2}$', fontsize=18)
    ax.set_ylabel('Counts', fontsize=18)

    if fax is None:
        plt.show()
    else:
        return


def ring_plot(data, fax=None, block=True):
    """Plots the Fabry Perot ring image
    
    Args:
        data (np.ndarray): 2D image data
        fax (tuple, optional): figure and axis
        block (bool): True if you want plt.show(block=True)
    """
    if fax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig, ax = fax

    cb = ax.imshow(data, cmap='Greys_r', origin='lower', interpolation=None)
    #cb = ax.imshow(data, cmap='Greys_r', interpolation=None)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cb, cax=cax, extend='max')
    if fax is not None:
        return
    else:
        plt.show(block=block)
