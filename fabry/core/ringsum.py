from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from numba import jit
'''
Core module contains ringsum codes that are the basis of this
analysis. 

Functions:
    get_binarr: returns binarr for ringsums
    smAng_ringsum: main ringsum function, uses small angle approx.
    proper_ringsum: additional ringsum function that makes no approx.
    locate_center: center finding function 
    new_ringsum: best ringsum to use currently
'''


def _histogram(bins, R, weights, out=None, label=None):
        '''
        Helper function to perform a histogram using multiprocessing

        Args:
            bins (np.array): bins in radius pixel space (does not include origin)
            R (np.ndarray): radius matrix for the camera
            weights (np.ndarray): pixel values for the camera
            out (mp.Queue): mulitprocessing.Queue to write data to
            label (str): label for the quadrant being ringsummed

        Returns:
            None if used with an output Queue
            sig (np.array): ringsum array
            label (str): label for the quadrant
        '''
        sig, _ = np.histogram(R, bins=np.concatenate((np.array([0.]),
                                bins)),weights=weights)
        if out and label:
            out.put((label, sig))
        else:
            return sig, label

def get_bin_edges(dat, x0, y0, binsize=0.1):
    '''
    returns equal area annuli bins for a ring image
    given a center and binsize

    Args:
        dat (np.ndarray): pixel values for the camera image
        x0 (float): x location (pixels) of the center of the image
        y0 (float): y location (pixels) of the center of the image
        binsize (float, default=0.1): smallest radial bin size 
                (occurs at the largest bin)

    Returns:
        binarr (np.array): bins for the ring sum in pixels
    '''
    # ny, nx = dat.shape
    # x = np.arange(0, nx, 1)
    # y = np.arange(0, ny, 1)

    # xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    # R = np.sqrt(xx**2 + yy**2)

    # xmin = xx[0, 0]
    # xmax = xx[0, -1]
    # ymin = yy[0, 0]
    # ymax = yy[-1, 0]

    # ri = int(np.min(np.abs([xmin, xmax, ymin, ymax])))
    # imax = int((ri ** 2. - 2. * ri - 1.) / (1. + 2. * ri) / binsize)

    # norm_radius = np.sqrt(2*binsize*ri + binsize**2)
    ny, nx = dat.shape
    x = np.array([0, nx-1]) - x0
    y = np.array([0, ny-1]) - y0
    xmin = np.abs(x).min()
    ymin = np.abs(y).min()
    ri = np.min([xmin, ymin])

    imax = int(np.floor(ri**2 / (2*ri-binsize) / binsize))
    i =  np.linspace(0, imax, imax+1)
    redges = np.sqrt(i * (2*ri - binsize)*binsize)

    return redges


def smAng_ringsum(dat, x0, y0, binsize=0.1, quadrants=False):
    '''
    Performs a constant area ring sum in pixel space. This is valid
    when the small angle approximation is taken for Cos(Theta). 

    Args:
        dat (np.ndarray): pixel values for the camera image
        x0 (float): x location (pixels) of the center of the image
        y0 (float): y location (pixels) of the center of the image
        binsize (float, default=0.1): smallest radial bin size 
                (occurs at the largest bin)
        quadrants (bool, default=False): if True, returns separate
                ringsums for each quadrant (used for center finding)

    Returns:
        binarr (np.array): bins for the ring sum in pixels
        ringsum (np.array): ringsum(s) If quadrants is True, there
                will be 4 separate ringsums in order: UL UR BL BR
    '''
    ny, nx = dat.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)

    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)

    xmin = xx[0, 0]
    xmax = xx[0, -1]
    ymin = yy[0, 0]
    ymax = yy[-1, 0]
    #print(xmin, xmax, ymin, ymax)
    ri = int(np.min(np.abs([xmin, xmax, ymin, ymax])))
    imax = int((ri ** 2. - 2. * ri - 1.) / (1. + 2. * ri) / binsize)

    norm_radius = np.sqrt(2*binsize*ri + binsize**2)
    binarr = np.sqrt(range(1, imax))*norm_radius
    
    xi0 = int(x0)
    yi0 = int(y0)

    i1 = [yi0-ri, yi0-ri, yi0, yi0]
    i2 = [yi0+1, yi0+1, yi0+ri+1, yi0+ri+1]
    j1 = [xi0-ri, xi0, xi0-ri, xi0]
    j2 = [xi0+1, xi0+ri+1, xi0+1, xi0+ri+1]

    procs = []
    nprocs = 4
    sigs = {}
    out = mp.Queue()
    labels = ['UL', 'UR', 'BL', 'BR']
    for k in range(nprocs):
        p = mp.Process(target=_histogram, args=(binarr, 
            R[i1[k]:i2[k], j1[k]:j2[k]], dat[i1[k]:i2[k], j1[k]:j2[k]]),
            kwargs={'out': out, 'label': labels[k]})
        procs.append(p)
        p.start()

    for i in range(nprocs):
        tup = out.get()
        sigs[tup[0]] = tup[1]

    for p in procs:
        p.join()

    # Return all 4 quadrant rings sums or sum them into one
    if quadrants:
        return binarr, sigs['UL'], sigs['UR'], sigs['BL'], sigs['BR']
    else:
        return binarr, sigs['UL']+sigs['UR']+sigs['BL']+sigs['BR']


def proper_ringsum(R, weights, m, L, d, peaks, lambda_min, lambda_max, delta_lambda, ndl=512):
    """
    Performs a proper ringsum with constant lambda spacing for each peak.
    This is quite slow and not implemented anywhere in the code currently,
    but has been included for completeness.

    Args:
        R (np.array): flattened array matrix from camera
        weights: flattened pixel value matrix from camera
        m (float): highest order number (not an integer value)
        L (float): camera focal length in pixels
        d (float): etalon spacing in mm
        peaks (list): location of peaks in pixels
        lambda_min (float): minimum wavelength in array
        lambda_max (float): maximum wavelength in array
        delta_lambda (float): delta wavelength spacing

    Returns:
        ringsum (list): a list of ringsums for each peak (order) in peaks
    """
    ringsum = []
    for idx, peak in enumerate(peaks):
        if isinstance(m, list):
            mm = m[idx]
        else:
            mm = m - idx
        rmin = np.sqrt((2 * d * L * 1.e6 / (mm * lambda_max)) ** 2 - L ** 2)
        rmax = np.sqrt((2 * d * L * 1.e6 / (mm * lambda_min)) ** 2 - L ** 2)
        inds = np.where(np.logical_and(R > rmin, R < rmax))[0]

        r0 = 2.*d*1.e6*L / mm / delta_lambda
        dr = np.sqrt((1. / (1. / np.sqrt(L ** 2 + rmin ** 2) - 1. / r0)) ** 2 - L ** 2) - rmin
        rR = [rmin, rmin + dr]
        j = 0
        while len(rR) <= ndl: #rR[-1] <= rmax:
            j += 1
            dr = np.sqrt((1. / (1. / np.sqrt(L ** 2 + rR[-1] ** 2) - 1. / r0)) ** 2 
                    - L ** 2) - rR[-1]
            rR += [rR[-1] + dr]

        sig, _ = np.histogram(R[inds], bins=rR, weights=weights[inds])
        # to correspond to my calibration wavelength array, sig is actually backwards...
        sig = sig[::-1]
        ringsum += [sig]
    return ringsum


def locate_center(data_in, xguess=None, yguess=None, maxiter=25, binsize=0.1, plotit=False, block_center=False):
    '''
    Finds the center of a ring pattern image by preforming ringsums.

    Args:
        data (np.ndarray): pixel values for the camera image
        xguess (float, default=None): guess of x location of center, if None, takes data center
        yguess (float, default=None): guess of y location of center, if None, takes data center
        maxiter (int, default=25): maximum number of center finding iterations
        binsize (float, default=0.1): smallest radial binsize (occurs at the largest bin)
        plotit (bool, default=False): plot fitting curves
        block_center (bool, default=False): block center of image (600x600 on xguess,yguess) when           finding center. This helps when there is a ring at the center of the image.

    Returns:
        x0 (float): x location of center
        y0 (float): y location of center
    '''
    if xguess is None:
        xguess = data.shape[1]/2.
    if yguess is None:
        yguess = data.shape[0]/2.
    print(xguess, yguess)    
    if block_center:
        data = np.copy(data_in)
        data[int(yguess-300):int(yguess+301),int(xguess-300):int(xguess+301)] = 0.0
    else:
        data = data_in

    if plotit:
        plt.ion()
        fig, ax = plt.subplots()
        line1 = None
        line2 = None
        line3 = None
        line4 = None

    print( "Center finding:")
    print( "start x0: {0} y0: {1}".format(xguess, yguess))
    for ii in range(maxiter):
        #print(xguess, yguess)
        #binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = smAng_ringsum(data, xguess, yguess, binsize=binsize,quadrants=True)
        binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = ringsum(data, xguess, yguess, binsize=binsize,quadrants=True)
        #jj = np.where(binarr < 300)
        #ULsigarr[jj] = 0.0
        #URsigarr[jj] = 0.0
        #BRsigarr[jj] = 0.0
        #BLsigarr[jj] = 0.0
        #fig, ax = plt.subplots()
        #ax.plot(binarr, ULsigarr)
        #ax.plot(binarr, URsigarr)
        #plt.show()
        thres = 0.3* np.max(ULsigarr + URsigarr)
        i = np.where(ULsigarr + URsigarr > thres)[0]

        # A cheat if the i goes to the edge of the array.  Makes sliding them past each other hard
        if len(ULsigarr) - i.max() < 60:
            i = i[0:-50]
        ni = len(i)
        #print(ni)
        ns = 25
        # ns = 13
        sarr = 2 * np.arange(-ns, ns+1, 1)
        sarr_max = np.max(sarr)
        UB = np.zeros(len(sarr))
        RL = np.zeros(len(sarr))

        for idx, ix in enumerate(sarr):
            UB[idx] = np.sum((ULsigarr[i-ix]+URsigarr[i-ix] - BLsigarr[i]-BRsigarr[i])**2) / (1.0*ni - np.abs(ix))
            RL[idx] = np.sum((URsigarr[i-ix]+BRsigarr[i-ix] - ULsigarr[i]-BLsigarr[i])**2) / (1.0*ni - np.abs(ix))
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
            RLcent = -RLfit[1] / (2*RLfit[0])

        # Dont jump fartther than max(sarr)
        if np.abs(RLcent) > sarr_max:
            RLcent = np.sign(RLcent) * np.max(sarr)

        if UBfit[0] < 0.0:
            # concave down
            UBcent = -2 * np.max(sarr) * np.sign(UBfit[1])
        else:
            # concave up
            UBcent = -UBfit[1] / (2*UBfit[0])

        # Dont jump fartther than max(sarr)
        if np.abs(UBcent) > sarr_max:
            UBcent = np.sign(UBcent) * np.max(sarr)

        if False:#~np.isfinite(RLcent):
            xguess -= binsize
        else:
            xguess -= RLcent * binsize

        if False:#~np.isfinite(UBcent):
            yguess += binsize
        else:
            yguess += UBcent * binsize
        
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
                leg = ax.legend()

            fig.canvas.draw()
            #plt.pause(0.0001)
            plt.pause(1.0)

        if np.sqrt((UBcent*binsize)**2 + (RLcent*binsize)**2) / binsize < 0.1:
            break
    if plotit:
        plt.close(fig)
        plt.ioff()
    return xguess, yguess


def new_ringsum(data, redges, x0, y0, use_weighted=False):
    """
    redges are all of the right edges (the first left edge is zero)
    """
    ny, nx = data.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)

    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)

    R = R.flatten()
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
        portion = slice(start,start+iedge)
        #if idx == 5000:
        #    fig, ax = plt.subplots()
        #    ax.hist(R[portion])
        #    plt.show()
        if use_weighted:
            means[idx], sigmas[idx] = calculate_weighted_mean(d[portion], np.sqrt(1.8*d[portion]))
        else:
            means[idx] = np.mean(d[portion])
            sigmas[idx] = np.std(d[portion]) / np.sqrt(len(d[portion]))

        lengths[idx] = len(d[portion])
        start += iedge
    #print(redges[-1] - redges[-2])
    #fig, ax = plt.subplots()
    #plt.plot(redges[0:-1]**2, np.diff(redges**2))
    #plt.show(block=False)

    #fig, ax = plt.subplots()
    #plt.plot(np.diff(redges))
    #plt.show()
    # fig, ax =plt.subplots()
    # ax.hist(lengths, bins='auto', density='True')
    # ax.set_xlabel('Number of Points in a Ring')
    # ax.set_title("{0:d} +/- {1:d}".format(int(np.mean(lengths)), int(np.std(lengths))))
    # plt.show()
    return means, sigmas


@jit(nopython=True)
def calculate_weighted_mean(data, error):
    idx = np.where(error > 0.0)
    err = error[idx]
    d = data[idx]

    weights = 1.0 / err**2

    denominator = np.nansum(weights)

    sigma = np.sqrt(1.0 / denominator)

    numerator = np.nansum(d * weights)

    mean = numerator / denominator

    return mean, sigma


@jit(nopython=True)
def super_pixelate(data, npix=2):
    n, m = data.shape

    n_new = n // npix
    m_new = m // npix

    d = np.zeros((n_new, m_new))

    for i in xrange(n_new):
        for j in xrange(m_new):
            n_idx = slice(i*npix, (i+1)*npix)
            m_idx = slice(j*npix, (j+1)*npix)
            # d[i, j] = np.mean(data[n_idx, m_idx])
            d[i, j] = np.sum(data[n_idx, m_idx])
    return d


def ringsum(data, x0, y0, binsize=0.1, quadrants=False, use_weighted=False):

    ny, nx = data.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)

    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)

    #xmin = xx[0, 0]
    #xmax = xx[0, -1]
    #ymin = yy[0, 0]
    #ymax = yy[-1, 0]

    #ri = int(np.min(np.abs([xmin, xmax, ymin, ymax])))
    #imax = int(np.floor(ri**2 / (2*ri-binsize) / binsize))
    #i =  np.linspace(0, imax, imax+1)
    #redges = np.sqrt(i * (2*ri - binsize)*binsize)

    redges = get_bin_edges(data, x0, y0, binsize=binsize)
    ri = int(redges[-1])

    #print(ri, test[-1])
    #fig, ax = plt.subplots()
    #ax.plot(redges)
    #ax.plot(test)
    #plt.show(block=True)

    #ri = int(redges[-1])+1
    #print(ri, int(redges[-1]))

    xi0 = int(x0)
    yi0 = int(y0)

    i1 = [yi0-ri, yi0-ri, yi0, yi0]
    i2 = [yi0+1, yi0+1, yi0+ri+1, yi0+ri+1]
    j1 = [xi0-ri, xi0, xi0-ri, xi0]
    j2 = [xi0+1, xi0+ri+1, xi0+1, xi0+ri+1]

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
                kwargs={'out': out, 'label': labels[k], 'use_weighted':False})
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
        portion = slice(start,start+iedge)
        if use_weighted:
            means[idx], sigmas[idx] = calculate_weighted_mean(d[portion], np.sqrt(1.8*d[portion]))
        else:
            means[idx] = np.mean(d[portion])
            sigmas[idx] = np.std(d[portion]) / np.sqrt(len(d[portion]))

        lengths[idx] = len(d[portion])
        start += iedge


    if out and label:
        out.put((label, means, sigmas))
    else:
        return means, sigmas


