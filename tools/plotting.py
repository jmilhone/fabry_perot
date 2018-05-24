import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

def tableau20_colors():
    '''
    good set of colors for plotting
    
    Returns:
        tableau20 (list): list of r,g,b color values
    '''
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
    '''
    custom histogram function for plotting multinest output

    Args:
        ax (matplotlib.axes): axes to plot with
        data (np.ndarray): array of values to histogram
        bins (int, optional): number of bins in histogram
    '''
    if bins is not None:
        hist, bins = np.histogram(data, density=True, bins=bins)
    else:
        hist, bins = np.histogram(data, density=True, bins='auto')

    bw = bins[1]-bins[0]
    
    if horizontal:
        ax.barh(bins[0:-1], hist*bw, height=bw)
        if data.max() > 1000:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        else:
            ax.get_yaxis().get_major_formatter().set_scientific(True)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
    else:
        ax.bar(bins[0:-1], hist*bw, width=bw)
        if data.max() > 1000:
            # I don't think this works
            #ax.get_xaxis().get_major_formatter().set_scientific(True)
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        else:
            ax.get_xaxis().get_major_formatter().set_scientific(True)

        ax.get_xaxis().get_major_formatter().set_useOffset(False)
    return bins

def my_hist2d(ax, data1, data2, bins=None, z=30):
    '''
    custom 2d histrogram option for plotting correlations

    Args:
        ax (matplotlib.axes): axes to plot with
        data1 (np.ndarray): array of values to histogram 
        data2 (np.ndarray): other array of values
        bins (int, optional): number of bins in histogram
    '''
    if bins is not None:
        hist, xx, yy = np.histogram2d(data1, data2, normed=True, bins=bins)
    else:
        hist, xx, yy = np.histogram2d(data1, data2, normed=True)

    dx = xx[1]-xx[0]
    dy = yy[1]-yy[0]
   
    im = ax.contourf(yy[0:-1], xx[0:-1], hist*dx*dy, z)
    #ax_divider = make_axes_locatable(ax)
    #cax = ax_divider.append_axes("right", size="7%", pad="2%")
    cb = plt.colorbar(im, ax=ax)

    if data1.max() > 1000:
        #ax.get_yaxis().get_major_formatter().set_scientific(True)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    else:
        ax.get_yaxis().get_major_formatter().set_scientific(False)
    if data2.max() > 1000:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #ax.get_xaxis().get_major_formatter().set_scientific(True)
    else:
        ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    return cb

class ClickBuilder():
    '''
    Little class for interactive plot features

    Returns x, y data from click(s)
    '''
    def __init__(self, fig):
        self.x = []
        self.y = []
        self.fig = fig
        self.cid = fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        self.x.append(event.xdata)
        self.y.append(event.ydata)

def center_plot(data, x0=None, y0=None):
    '''
    Runs clickable center plot for an image

    Args:
        data (np.ndarray): pixel values for image
        x0 (float, default=None): x center for zoomed image, if
            none will use the center of data array
        y0 (float, default=None): y center for zoomed image, if
            none will use the cneter of data array
    Returns:
        x (float): x value of click guess
        y (float): y value of click guess
    '''
    if x0 is None:
        x0 = data.shape[1]/2.
    if y0 is None:
        y0 = data.shape[0]/2.
    
    fig, axs = plt.subplots(1,2,figsize=(12,6))
    cb = axs[0].imshow(data, cmap='gray', origin='lower')
    axs[0].add_patch(Rectangle((x0-500,y0-500),1001,1001,lw=2,linestyle='--',ec='red',fc='none'))
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(cb, cax=cax)
    axs[1].imshow(data,cmap='gray',origin='lower')
    axs[1].set_title('Click for center guess')
    plt.setp(axs[1].spines.values(), color='red')
    [i.set_linewidth(5) for i in axs[1].spines.itervalues()]
    [i.set_linestyle('--') for i in axs[1].spines.itervalues()]
    axs[1].set_ylim(y0-500,y0+501)
    axs[1].set_xlim(x0-500,x0+501)

    center_guess = ClickBuilder(fig)
    plt.tight_layout()
    plt.show()
    fig.canvas.mpl_disconnect(center_guess.cid)
    return center_guess.x[0], center_guess.y[0]
    
def ringsum_click(r, sig, title='Click Me!'):
    '''
    Runs clickable ringsum plot

    Args:
        r (np.ndarray): bin values
        sig (np.ndarray): ringsum values
        title (str): instructions for clicker displayed as plot title

    Returns:
        x (list): x value(s) of click(s)
        y (list): y value(s) of click(s)
    '''
    fig, axs = plt.subplots(figsize=(10,6))
    axs.plot(r, sig, lw=2)
    axs.set_title(title, fontsize=22)
    axs.set_xlabel('R (pixels)', fontsize=18)
    axs.set_ylabel('Counts', fontsize=18)
    clicks = ClickBuilder(fig)
    plt.show()
    fig.canvas.mpl_disconnect(clicks.cid)
    return clicks.x, clicks.y

def peak_plot(r, sig, peaks, orders, fax=None, pkerr=None):
    '''
    plots ringsum data with labeled peaks and orders

    Args:
        r (np.ndarray): r array from ringsum
        sig (np.ndarray): ringsum signal
        peaks (dict): dictionary of peak locations
            the keys are the wavelengths
        orders (dict): dictionary of peak orders
            the keys are the wavelengths
        fax (tuple, default=None): the figure and axis
            handles for the plot, if adding to existing
            plot. Default is None, which will make a 
            new set of figure and axes and run plt.show()
        pkerr (dict, default=None): optional error dictionary
            which will plot the axvspan of the peak error
    '''
    if fax is None:
        fig, ax = plt.subplots(figsize=(10,6))
    else:
        fig, ax = fax

    colors = tableau20_colors()
    ax.plot(r, sig, color=colors[0], lw=2)
    i = 1
    for key in peaks.keys():
        for j, pk in enumerate(peaks[key]):
            if pkerr is None:
                ax.axvline(pk, lw=2, color=colors[i], label='{0}: j={1}'.format(key, orders[key][j]))
            else:
                ax.axvspan(pk+pkerr[key][j], pk-pkerr[key][j], color=colors[i], label='{0}: j={1}'.format(key, orders[key][j]), alpha=0.7)
            i += 1
            if i > len(colors):
                i = 0
    ax.legend(fontsize=14)
    ax.set_xlabel('R (pixels)',fontsize=18)
    ax.set_ylabel('Counts', fontsize=18)
    
    if fax is None:
        plt.show()
    else:
        return

def ring_plot(data, fax=None, block=True):
    if fax is None:
        fig, ax = plt.subplots(figsize=(10,8))
    else:
        fig, ax = fax
    
    cb = ax.imshow(data, cmap='Greys_r', origin='lower', interpolation=None)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cb, cax=cax, extend='max')
    if fax is not None:
        return
    else:
        plt.show(block=block)
