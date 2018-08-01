from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from numba import jit

"""
Core module contains ringsum codes that are the basis of this
analysis.

Functions:
    get_binarr: returns binarr for ringsums
    smAng_ringsum: main ringsum function, uses small angle approx.
    proper_ringsum: additional ringsum function that makes no approx.
    locate_center: center finding function 
    new_ringsum: best ringsum to use currently
"""


def get_bin_edges(dat, x0, y0, binsize=0.1):
    """Returns equal area annuli bins for a ring image given a center and binsize
    
    .. math::
        r_n = \sqrt{n \left(2 r_N - \Delta  \\right) \Delta } \quad n \in 0,...,N
    
    .. math:: 
        \Delta = \\text{ binsize}

    Args:
        dat (np.ndarray): pixel values for the camera image
        x0 (float): x location (pixels) of the center of the image
        y0 (float): y location (pixels) of the center of the image
        binsize (float, optional): smallest radial bin size (occurs at the largest bin)

    Returns:
        np.ndarray: bins for the ring sum in pixels
    """

    ny, nx = dat.shape
    x = np.array([0, nx - 1]) - x0
    y = np.array([0, ny - 1]) - y0
    xmin = np.abs(x).min()
    ymin = np.abs(y).min()
    ri = np.min([xmin, ymin])

    imax = int(np.floor(ri ** 2 / (2 * ri - binsize) / binsize))
    i = np.linspace(0, imax, imax + 1)
    redges = np.sqrt(i * (2 * ri - binsize) * binsize)

    return redges


def locate_center(data_in, xguess=None, yguess=None, maxiter=25, binsize=0.1, plotit=False, block_center=False,
                  printit=True):
    """
    Finds the center of a ring pattern image by preforming ringsums.

    Args:
        data_in (np.ndarray): pixel values for the camera image
        xguess (float): guess of x location of center, if None, takes data center
        yguess (float): guess of y location of center, if None, takes data center
        maxiter (int): maximum number of center finding iterations
        binsize (float): smallest radial binsize (occurs at the largest bin)
        plotit (bool): plot fitting curves
        block_center (bool): block center of image (600x600 on xguess,yguess) when finding center. This helps when
            there is a ring at the center of the image.
        printit (bool): print center find progress if True

    Returns:
        tuple (float, float): x and y location of the center
    """
    if xguess is None:
        xguess = data_in.shape[1] / 2.

    if yguess is None:
        yguess = data_in.shape[0] / 2.

    if printit:
        print(xguess, yguess)

    if block_center:
        data = np.copy(data_in)
        data[int(yguess - 300):int(yguess + 301), int(xguess - 300):int(xguess + 301)] = 0.0
    else:
        data = data_in

    line1 = None
    line2 = None
    fig = None
    ax = None

    if plotit:
        plt.ion()
        fig, ax = plt.subplots()

    if printit:
        print("Center finding:")
        print("start x0: {0} y0: {1}".format(xguess, yguess))

    for ii in range(maxiter):
        binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = ringsum(data, xguess, yguess, binsize=binsize, quadrants=True)

        thres = 0.3 * np.max(ULsigarr + URsigarr)
        i = np.where(ULsigarr + URsigarr > thres)[0]

        # A cheat if the i goes to the edge of the array.  Makes sliding them past each other hard
        if len(ULsigarr) - i.max() < 60:
            i = i[0:-50]
        ni = len(i)

        ns = 25
        # ns = 13
        sarr = 2 * np.arange(-ns, ns + 1, 1)
        sarr_max = np.max(sarr)
        UB = np.zeros(len(sarr))
        RL = np.zeros(len(sarr))

        for idx, ix in enumerate(sarr):
            UB[idx] = np.sum((ULsigarr[i - ix] + URsigarr[i - ix] - BLsigarr[i] - BRsigarr[i]) ** 2) / (
                    1.0 * ni - np.abs(ix))
            RL[idx] = np.sum((URsigarr[i - ix] + BRsigarr[i - ix] - ULsigarr[i] - BLsigarr[i]) ** 2) / (
                    1.0 * ni - np.abs(ix))
        RLfit = np.polyfit(sarr, RL, 2)
        UBfit = np.polyfit(sarr, UB, 2)

        """
        The logic for the center jump is matching A(x-x0)^2 = C0 x^2 + C1 x + C2
        x0 = - C1 / (2 C0)
        """
        if RLfit[0] < 0.0:
            # Concave down fit
            RLcent = -2 * np.max(sarr) * np.sign(RLfit[1])
        else:
            # concave up
            RLcent = -RLfit[1] / (2 * RLfit[0])

        # Dont jump fartther than max(sarr)
        if np.abs(RLcent) > sarr_max:
            RLcent = np.sign(RLcent) * np.max(sarr)

        if UBfit[0] < 0.0:
            # concave down
            UBcent = -2 * np.max(sarr) * np.sign(UBfit[1])
        else:
            # concave up
            UBcent = -UBfit[1] / (2 * UBfit[0])

        # Dont jump fartther than max(sarr)
        if np.abs(UBcent) > sarr_max:
            UBcent = np.sign(UBcent) * np.max(sarr)

        if False:  # ~np.isfinite(RLcent):
            xguess -= binsize
        else:
            xguess -= RLcent * binsize

        if False:  # ~np.isfinite(UBcent):
            yguess += binsize
        else:
            yguess += UBcent * binsize

        if printit:
            print("{2:d}, update x0: {0}, y0: {1}".format(xguess, yguess, ii))

        if plotit:
            if line1 is not None:
                line1.set_data(sarr, UB)
            else:
                line1, = ax.plot(sarr, UB, 'r', lw=1, label='UD')

            if line2 is not None:
                line2.set_data(sarr, RL)
            else:
                line2, = ax.plot(sarr, RL, 'b', lw=1, label='RL')

            if ii == 0:
                ax.legend()

            fig.canvas.draw()
            plt.pause(1.0)

        if np.sqrt((UBcent * binsize) ** 2 + (RLcent * binsize) ** 2) / binsize < 0.1:
            break

    if plotit:
        plt.close(fig)
        plt.ioff()

    return xguess, yguess


@jit(nopython=True)
def calculate_weighted_mean(data, error):
    """Calculates the weighted mean of data with standard deviation error

    Args:
        data (np.ndarray): data array
        error (np.ndarray): standard deviation error bar for data

    Returns:
        tuple (float, float): weighted mean and weighted standard deviation
    """
    idx = np.where(error > 0.0)
    err = error[idx]
    d = data[idx]

    weights = 1.0 / err ** 2

    denominator = np.nansum(weights)

    sigma = np.sqrt(1.0 / denominator)

    numerator = np.nansum(d * weights)

    mean = numerator / denominator

    return mean, sigma


@jit(nopython=True)
def super_pixelate(data, npix=2):
    """Creates super pixels for image data

    Args:
        data (np.ndarray): 2d image data

        npix (int, optional): integer number of pixels to create an npix x npix super pixel

    Returns:
        np.ndarray: New image made from the super pixels
    """
    n, m = data.shape

    n_new = n // npix
    m_new = m // npix

    d = np.zeros((n_new, m_new))

    for i in xrange(n_new):
        for j in xrange(m_new):
            n_idx = slice(i * npix, (i + 1) * npix)
            m_idx = slice(j * npix, (j + 1) * npix)
            # d[i, j] = np.mean(data[n_idx, m_idx])
            d[i, j] = np.sum(data[n_idx, m_idx])
    return d


def ringsum(data, x0, y0, binsize=0.1, quadrants=False, use_weighted=False):
    """Returns a equal annulus area ringsum centered at (x0, y0) from data

    Args:
        data (np.ndarray): 2d image data
        x0 (float): center location in x
        y0 (float): center location in y
        binsize (float, optional): the delta r of the last annulus, default=0.1
        quadrants (bool): split the ringsum into 4 quadrants to use multiprocessing, default=False
        use_weighted (bool): use a weighted mean, default=False

    Returns:
        tuple
    """
    ny, nx = data.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)

    xx, yy = np.meshgrid(1. * x - x0, 1. * y - y0)
    R = np.sqrt(xx ** 2 + yy ** 2)

    redges = get_bin_edges(data, x0, y0, binsize=binsize)
    ri = int(redges[-1])

    xi0 = int(x0)
    yi0 = int(y0)

    i1 = [yi0 - ri, yi0 - ri, yi0, yi0]
    i2 = [yi0 + 1, yi0 + 1, yi0 + ri + 1, yi0 + ri + 1]
    j1 = [xi0 - ri, xi0, xi0 - ri, xi0]
    j2 = [xi0 + 1, xi0 + ri + 1, xi0 + 1, xi0 + ri + 1]

    rarr = 0.5 * (redges[0:-1] + redges[1:])

    if quadrants:
        procs = []
        nprocs = 4
        sigs = {}
        out = mp.Queue()
        labels = ['UL', 'UR', 'BL', 'BR']
        for k in range(nprocs):
            p = mp.Process(target=_ringsum, args=(redges[1:],
                                                  R[i1[k]:i2[k], j1[k]:j2[k]], data[i1[k]:i2[k], j1[k]:j2[k]]),
                           kwargs={'out': out, 'label': labels[k], 'use_weighted': False})
            procs.append(p)
            p.start()

        for i in range(nprocs):
            tup = out.get()
            sigs[tup[0]] = tup[1]

        for p in procs:
            p.join()

        return rarr, sigs['UL'], sigs['UR'], sigs['BL'], sigs['BR']
    else:
        sig, sigma = _ringsum(redges[1:], R, data, out=None, label=None, use_weighted=False)
        return rarr, sig, sigma


def _ringsum(redges, radii, data, out=None, label=None, use_weighted=False):
    """Helper function for ringsumming

    Args:
        redges (np.ndarray): bin edges (does not include origin)
        radii (np.ndarray): radii to bin
        data (np.ndarray): image weights for the radii to be binned with
        out (mp.Queue, optional): multiprocessing queue to place results in if needed
        label (list, optional): label to put with results when placing in out
        use_weighted (bool): use weighted mean if True, default=False

    Returns:
        None if use_weighted
        tuple (np.ndarray, np.ndarray): ring sum, ring sum standard deviations
    """
    # redges does not include zero!
    n = len(redges)

    R = radii.flatten()
    d = data.flatten()

    indsort = np.argsort(R)
    R = R[indsort]
    d = d[indsort]

    n = len(redges)
    means = np.zeros(n)
    sigmas = np.zeros(n)
    lengths = np.zeros(n)
    start = 0
    for idx, edge in enumerate(redges):
        iedge = np.searchsorted(R[start:], edge, side='right')
        portion = slice(start, start + iedge)
        if use_weighted:
            means[idx], sigmas[idx] = calculate_weighted_mean(d[portion], np.sqrt(1.8 * d[portion]))
        else:
            means[idx] = np.mean(d[portion])
            sigmas[idx] = np.std(d[portion]) / np.sqrt(len(d[portion]))

        lengths[idx] = len(d[portion])
        start += iedge

    if out and label:
        out.put((label, means, sigmas))
    else:
        return means, sigmas
