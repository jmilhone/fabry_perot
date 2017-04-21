import numpy as np
import matplotlib.pyplot as plt
from image_helpers import get_image_data, quick_plot
from os.path import join
import fitting
import time
import multiprocessing as mp



def ringsum(R, weights, m0, L, d, peaks, delta_lambda, ndl=512):
    # L = np.float64(L)
    # d = np.float64(d)
    # peaks = np.array(peaks, dtype=np.float64)
    # delta_lambda = np.float64(delta_lambda)
    d_lam_arr = np.arange(-ndl/2, ndl/2+1, 1) * delta_lambda
    # print d_lam_arr
    ringsums = []
    for j, peak in enumerate(peaks):


        sq = np.sqrt(1.0 + (peak/L)**2)
        num = 2 * d * 1.e6 * sq
        denom = 2*d*1.e6 + d_lam_arr[::-1] * (m0-j) * sq * -1.0  # added -1 to fix signs
        temp = num / denom
        # print temp
        rarr = L * np.sqrt(temp**2 - 1.0)
        # print m0, rarr.min(), rarr.max(), d_lam_arr.max() - d_lam_arr.min()
        # rarr = rarr[::-1]
        rmin = np.nanmin(rarr)
        rmax = np.nanmax(rarr)
        # inds = np.where( np.logical_and(R >= rmin, R < rmax))

        # rmin = 600.0
        # rmax = 750.0
        # rarr = np.linspace(rmin, rmax, 100)
        inds = np.where( np.logical_and(R >= rmin, R < rmax))
        # print rmin, rmax, peak
        # print inds
        # print rmin, rmax
        # print len(inds)
        sig, _ = np.histogram(R[inds], bins=rarr[::-1], weights=weights[inds])
        # print len(sig), len(rarr), len(d_lam_arr)
        # plt.plot(rarr[0:-1], sig)
        # plt.plot(d_lam_arr[0:-1], sig[::-1], '.')
        # plt.plot(d_lam_arr, rarr)
        # plt.show()
        # ringsums += [y(sig[::-1])]

        ringsums.append(sig[::-1])
        # plt.plot(R[inds], weights[inds])
        # print rarr
        # plt.plot(np.diff(rarr))


        # fig, ax = plt.subplots()
        # ax.imshow(weights)
        # ax.set_aspect(1.0)
        # theta = np.linspace(0, 2*np.pi, 1000)
        # plt.plot(rmin*np.cos(theta)+3018.5, rmin*np.sin(theta)+2010.5, 'r')
        # plt.plot(rmax*np.cos(theta)+3018.5, rmax*np.sin(theta)+2010.5, 'g')
        # plt.show()
        # print rarr
    print type(ringsums), type(ringsums[0])
    return ringsums, d_lam_arr[0:-1]

def proper_ringsum(R, weights, m, L, d, peaks, lambda_min, lambda_max, delta_lambda, ndl=512):
    """
    Performs a proper ringsum with constant lambda spacing for each peak.

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
            dr = np.sqrt((1. / (1. / np.sqrt(L ** 2 + rR[-1] ** 2) - 1. / r0)) ** 2 - L ** 2) - rR[-1]
            # print j, dr, rR[-1], rmax
            rR += [rR[-1] + dr]
        # print idx
        # print len(rR)
        # print rmin, rR[0]
        # print rmax, rR[-1], "\n"

        sig, _ = np.histogram(R[inds], bins=rR, weights=weights[inds])
        # to correspond to my calibration wavelength array, sig is actually backwards...
        sig = sig[::-1]
        ringsum += [sig]
    return ringsum
    # ny, nx = data.shape
    # x = np.arange(1, nx+1, 1)
    # y = np.arange(1, ny+1, 1)
    #
    # xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    # R = np.sqrt(xx**2 + yy**2)


def _histogram(bins, R, weights, out=None, label=None):

    """
    Helper function for quick_ringsum to perform a histogram
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
    """
    sig, _ = np.histogram(R, bins=np.concatenate((np.array([0.]), bins)), weights=weights)
    if out and label:
        out.put((label, sig))
    else:
        return sig, label

def quick_ringsum(dat, x0, y0, binsize=0.1, quadrants=True):

    """
    Performs a constant area ring sum in pixel space.  Returns ring sums for
     the 4 qudrants if quadrants is True

    Args:
        dat (np.ndarray): pixel values for the camera image
        x0 (float): x location of the center of image
        y0 (float): y location of the center of image
        binsize (float): smallest radial bin size (occurs at the largest bin)
        quadrants (bool): True if quadrant ring sums are returned

    Returns:
        binarr (np.array): bins for the ring sum in pixels
        ringsum (np.array): ring sums (there are 4 of them if quadrants=True,
         UL UR BL Br for the order)
    """
    ny, nx = dat.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)

    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)

    xmin = xx[0, 0]
    xmax = xx[0, -1]
    ymin = yy[0, 0]
    ymax = yy[-1, 0]

    ri = int(np.min(np.abs([xmin, xmax, ymin, ymax])))
    imax = int((ri ** 2. - 2. * ri - 1.) / (1. + 2. * ri) / binsize)

    # """Ken's way of building the bins"""
    # binarr = np.fromfunction(lambda i: np.sqrt(2. * (i + 1.) * ri * binsize + (i + 1.) * binsize ** 2.), (imax,),
    #                          dtype='float64')
    norm_radius = np.sqrt(2*binsize*ri + binsize**2)
    binarr = np.sqrt(range(1, imax))*norm_radius

    # xi0 = int(round(x0))
    # yi0 = int(round(y0))

    xi0 = int(x0)
    yi0 = int(y0)

    # i1 = [0, 0, yi0, yi0]
    # i2 = [yi0+1, yi0+1, ny, ny]
    # j1 = [0, xi0, 0, xi0]
    # j2 = [xi0+1, nx, xi0+1, nx]

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
        # print R[i1[k]:i2[k], j1[k]:j2[k]].shape
        p = mp.Process(target=_histogram, args=(binarr, R[i1[k]:i2[k], j1[k]:j2[k]], dat[i1[k]:i2[k], j1[k]:j2[k]]),
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


def locate_center(data, xguess, yguess, maxiter=25, binsize=0.1, plotit=False):
    """
    Finds the center of the image by performing rings.

    Args:
        data (np.ndarray): pixel values for the camera image
        xguess (float): guess of x location of the center of image
        yguess (float): guess of y location of the center of image
        maxiter (int): maximum number of iterations for center finding
        binsize (float): smallest radial bin size (occurs at the largest bin)
        plotit (bool): plots minimization for each iteration if True

    Returns:
        x0 (float): x location of center
        y0 (float): y location of center
    """
    if plotit:
        plt.ion()
        fig, ax = plt.subplots()
        line1 = None
        line2 = None
        line3 = None
        line4 = None

    print "Center finding:"
    print "start x0: {0} y0: {1}".format(xguess, yguess)
    for ii in range(maxiter):

        # t0 = time.time()
        binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = quick_ringsum(data, xguess, yguess, binsize=binsize)
        # print time.time()-t0

        thres = 0.3* np.max(ULsigarr + URsigarr)
        i = np.where(ULsigarr + URsigarr > thres)[0]

        # A cheat if the i goes to the edge of the array.  Makes sliding them past each other hard
        if len(ULsigarr) - i.max() < 60:
            i = i[0:-50]
        ni = len(i)
        #ns = 25
        ns = 11
        sarr = 2 * np.arange(-ns, ns+1, 1)
        sarr_max = np.max(sarr)
        UB = np.zeros(len(sarr))
        RL = np.zeros(len(sarr))

        # sarr = np.arange(-2*ns, 2*ns + 1, 1)
        for idx, ix in enumerate(sarr):
            # print i-ix, i
            UB[idx] = np.sum((ULsigarr[i-ix]+URsigarr[i-ix] - BLsigarr[i]-BRsigarr[i])**2) / (1.0*ni - np.abs(ix))
            RL[idx] = np.sum((URsigarr[i-ix]+BRsigarr[i-ix] - ULsigarr[i]-BLsigarr[i])**2) / (1.0*ni - np.abs(ix))
            #UB[idx] = np.sum((ULsigarr[i-ix]+URsigarr[i-ix] - BLsigarr[i+ix]-BRsigarr[i+ix])**2) / (1.*ni-np.abs(ix))
        RLfit = np.polyfit(sarr, RL, 2)
        UBfit = np.polyfit(sarr, UB, 2)


        """
        The logic for the center jump is matching A(x-x0)^2 = C0 x^2 + C1 x + C2
        x0 = - C1 / (2 C0)
        """
        if RLfit[0] < 0.0:
            # print "RL concave down"
            # Concave down fit
            RLcent = -2 * np.max(sarr) * np.sign(RLfit[1])
        else:
            # concave up
            RLcent = -RLfit[1] / (2*RLfit[0])

        # Dont jump fartther than max(sarr)
        if RLcent > sarr_max:
            RLcent = np.sign(RLcent) * np.max(sarr)


        if UBfit[0] < 0.0:
            # concave down
            # print "UB concave down"
            UBcent = -2 * np.max(sarr) * np.sign(UBfit[1])
        else:
            # concave up
            UBcent = -UBfit[1] / (2*UBfit[0])

        # Dont jump fartther than max(sarr)
        if RLcent > sarr_max:
            UBcent = np.sign(UBcent) * np.max(sarr)

        # Rough conversion to pixels
        # yguess -= UBcent * binsize
        yguess += UBcent * binsize
        xguess -= RLcent * binsize

        # yguess -= 1
        # print "subtracting 1 from yguess for no reason other than trying to match Cooper"

        print "{2:d}, update x0: {0}, y0: {1}".format(xguess, yguess, ii)
        print "{0}, {1}".format(UBcent, RLcent)
        if plotit:
            # fig, ax = plt.subplots()
            if line1 is not None:
                line1.set_data(sarr, UB)
            else:
                line1, = ax.plot(sarr, UB, 'r', lw=1, label='UD')

            if line2 is not None:
                line2.set_data(sarr, RL)
            else:
                line2, = ax.plot(sarr, RL, 'b', lw=1, label='RL')

            # if line3 is not None:
            #     line3.set_data(sarr, np.polyval(UBfit, sarr))
            # else:
            #     line3, = ax.plot(sarr, np.polyval(UBfit, sarr), '--r', lw=1, label='fit UD')
            #
            # if line4 is not None:
            #     line4.set_data(sarr, np.polyval(RLfit, sarr))
            # else:
            #     line4, = ax.plot(sarr, np.polyval(RLfit, sarr), '--b', lw=1, label='fit RL')
            if ii == 0:
                leg = ax.legend()

            fig.canvas.draw()
            plt.pause(0.0001)
            # plt.show()

        if np.sqrt((UBcent*binsize)**2 + (RLcent*binsize)**2) / binsize < 0.1:
            break
        # U = ULsigarr + URsigarr
        # B = BLsigarr + BRsigarr
        # R = URsigarr + BRsigarr
        # L = ULsigarr + BLsigarr
        # plt.plot(binarr, U)
        # plt.show()
        # Upeaks = fitting.peak_and_fit(binarr, U, thres=.6, npeaks=1)
        # Bpeaks = fitting.peak_and_fit(binarr, B, thres=.6, npeaks=1)
        # Rpeaks = fitting.peak_and_fit(binarr, R, thres=.6, npeaks=1)
        # Lpeaks = fitting.peak_and_fit(binarr, L, thres=.6, npeaks=1)
        #
        # UDloc = 0.5 * (Upeaks[0] + Bpeaks[0])
        # RLloc = 0.5 * (Rpeaks[0] + Lpeaks[0])
        #
        # print xguess, RLloc, 0.5*(RLloc - Rpeaks[0])
        # print yguess, UDloc, 0.5*(UDloc - Upeaks[0])
        #
        # xguess -= 0.5*(RLloc - Upeaks[0])
        # yguess += 0.5*(UDloc - Rpeaks[0])
        #
        # fig, ax = plt.subplots(2)
        # ax[0].plot(binarr, U)
        # ax[0].plot(binarr, B)
        # ax[1].plot(binarr, R)
        # ax[1].plot(binarr, L)
        # plt.show()


        # for idx, ix in enumerate(sarr):
        #     UB[idx] = np.sum((ULsigarr[i-ix]+URsigarr[i-ix] - BLsigarr[i+ix]-BRsigarr[i+ix])**2) / (1.*ni-np.abs(ix))
        #     RL[idx] = np.sum((URsigarr[i-ix]+BRsigarr[i-ix] - ULsigarr[i+ix]-BLsigarr[i+ix])**2) / (1.*ni-np.abs(ix))

        # RLfit = np.polyfit(sarr, RL, 2)
        # UBfit = np.polyfit(sarr, UB, 2)


        # """
        # The logic for the center jump is matching A(x-x0)^2 = C0 x^2 + C1 x + C2
        # x0 = - C1 / (2 C0)
        # """
        # if RLfit[0] < 0.0:
        #     # print "RL concave down"
        #     # Concave down fit
        #     RLcent = -2 * np.max(sarr) * np.sign(RLfit[1])
        # else:
        #     # concave up
        #     RLcent = -RLfit[1] / (2*RLfit[0])

        # # Dont jump fartther than max(sarr)
        # if RLcent > sarr_max:
        #     RLcent = np.sign(RLcent) * np.max(sarr)


        # if UBfit[0] < 0.0:
        #     # concave down
        #     # print "UB concave down"
        #     UBcent = -2 * np.max(sarr) * np.sign(UBfit[1])
        # else:
        #     # concave up
        #     UBcent = -UBfit[1] / (2*UBfit[0])

        # # Dont jump fartther than max(sarr)
        # if RLcent > sarr_max:
        #     UBcent = np.sign(UBcent) * np.max(sarr)

        # # Rough conversion to pixels
        # # yguess -= UBcent * binsize
        # yguess += UBcent * binsize
        # xguess -= RLcent * binsize

#         # yguess -= 1
#         # print "subtracting 1 from yguess for no reason other than trying to match Cooper"

        # print "{2:d}, update x0: {0}, y0: {1}".format(xguess, yguess, ii)

        # if plotit:
        #     # fig, ax = plt.subplots()
        #     if line1 is not None:
        #         line1.set_data(sarr, UB)
        #     else:
        #         line1, = ax.plot(sarr, UB, 'r', lw=1, label='UD')

            # if line2 is not None:
            #     line2.set_data(sarr, RL)
            # else:
            #     line2, = ax.plot(sarr, RL, 'b', lw=1, label='RL')

            # # if line3 is not None:
            # #     line3.set_data(sarr, np.polyval(UBfit, sarr))
            # # else:
            # #     line3, = ax.plot(sarr, np.polyval(UBfit, sarr), '--r', lw=1, label='fit UD')
            # #
            # # if line4 is not None:
            # #     line4.set_data(sarr, np.polyval(RLfit, sarr))
            # # else:
            # #     line4, = ax.plot(sarr, np.polyval(RLfit, sarr), '--b', lw=1, label='fit RL')
            # if ii == 0:
            #     leg = ax.legend()

            # fig.canvas.draw()
            # plt.pause(0.0001)
            # # plt.show()

        # if np.sqrt((UBcent*binsize)**2 + (RLcent*binsize)**2) / binsize < 0.1:
        #     break
    if plotit:
        plt.close(fig)
        plt.ioff()
    return xguess, yguess


if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    # Normal Ar plasma shot
    # shot_number = 15676
    # fname = join(folder, "{0:07d}_000.nef".format(shot_number))
    # bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))

    # Thorium Calibration with Ar 488 nm filter
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    data = get_image_data(fname, bg_fname, color='b')

    # quick_plot(data)
    ny, nx = data.shape
    # No initial guess
    x0 = nx/2.
    y0 = ny/2.
    # Very close guess to calibration center
    # x0, y0 = (3069.688, 2032.854)
    # plt.plot(data[y0, 0:x0])
    # plt.show()
    t0 = time.time()
    x0, y0 = locate_center(data, x0, y0, plotit=True)
    t1 = time.time()
    print t1-t0
    binarr, sigarr = quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    plt.plot(binarr, sigarr)
    plt.show()
    print time.time()-t1
    # plt.plot(binarr, sigarr)
    # plt.show()
    peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=True)



