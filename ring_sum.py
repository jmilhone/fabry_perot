import numpy as np
import matplotlib.pyplot as plt
from image_helpers import get_image_data
from os.path import join


def quick_ringsum(dat, x0, y0, binsize=0.1, quadrants=True):

    ny, nx = dat.shape
    x = np.arange(1, nx+1, 1)
    y = np.arange(1, ny+1, 1)

    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)

    xmin = xx[0, 0]
    xmax = xx[0, -1]
    ymin = yy[0, 0]
    ymax = yy[-1, 0]

    ri = int(np.min(np.abs([xmin, xmax, ymin, ymax])))
    imax = int((ri ** 2. - 2. * ri - 1.) / (1. + 2. * ri) / binsize)

    """Ken's way of building the bins"""
    # binarr = np.fromfunction(lambda i: np.sqrt(2. * (i + 1.) * ri * binsize + (i + 1.) * binsize ** 2.), (imax,),
    #                          dtype='float64')

    """This is Cooper's way of building the bins"""
    r = 0
    r_arr = []
    delta_r = [np.sqrt(2*binsize*ri + binsize**2)]
    while r < ri:
        # print r
        r += delta_r[-1]
        r_arr += [r]
        delta_r += [np.sqrt(r_arr[-1]**2 + 2*binsize*ri + binsize**2) - r_arr[-1]]
    binarr = np.array(r_arr)

    xi0 = int(round(x0))
    yi0 = int(round(y0))

    if quadrants:
        ULdata = dat[0:yi0 + 1, 0:xi0 + 1]
        URdata = dat[0:yi0+1, xi0:]
        BLdata = dat[yi0:, 0:xi0 + 1]
        BRdata = dat[yi0:, xi0:]

        RUL = R[0:yi0+1, 0:xi0+1]
        RUR = R[0:yi0+1, xi0:]
        RBL = R[yi0:, 0:xi0+1]
        RBR = R[yi0:, xi0:]

        ULsigarr, _ = np.histogram(RUL, bins=np.concatenate((np.array([0.]), binarr)), weights=ULdata)
        URsigarr, _ = np.histogram(RUR, bins=np.concatenate((np.array([0.]), binarr)), weights=URdata)
        BLsigarr, _ = np.histogram(RBL, bins=np.concatenate((np.array([0.]), binarr)), weights=BLdata)
        BRsigarr, _ = np.histogram(RBR, bins=np.concatenate((np.array([0.]), binarr)), weights=BRdata)

        return binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr

    else:
        sigarr, _ = np.histogram(R, bins=np.concatenate((np.array([0.]), binarr)), weights=dat)
        return binarr, sigarr


def locate_center(data, xguess, yguess, maxiter=25, binsize=0.1, plotit=False):
    for ii in range(maxiter):
        binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = quick_ringsum(data, xguess, yguess, binsize=binsize)

        thres = 0.35 * np.max(ULsigarr + URsigarr)
        i = np.where(ULsigarr + URsigarr > thres)[0]

        # A cheat if the i goes to the edge of the array.  Makes sliding them past each other hard
        if len(ULsigarr) - i.max() < 50:
            i = i[0:-50]
        ni = len(i)
        ns = 25
        sarr = 2 * np.arange(-ns, ns+1, 1)
        UB = np.zeros(len(sarr))
        RL = np.zeros(len(sarr))

        for idx, ix in enumerate(sarr):
            UB[idx] = np.sum((ULsigarr[i-ix]+URsigarr[i-ix] - BLsigarr[i+ix]-BRsigarr[i+ix])**2) / (1.*ni-np.abs(ix))
            RL[idx] = np.sum((URsigarr[i-ix]+BRsigarr[i-ix] - ULsigarr[i+ix]-BLsigarr[i+ix])**2) / (1.*ni-np.abs(ix))

        RLfit = np.polyfit(sarr, RL, 2)
        UBfit = np.polyfit(sarr, UB, 2)


        """
        The logic for the center jump is matching A(x-x0)^2 = C0 x^2 + C1 x + C2
        x0 = - C1 / (2 C0)
        """
        if RLfit[0] < 0.0:
            print "RL concave down"
            # Concave down fit
            RLcent = -2 * np.max(sarr) * np.sign(RLfit[1])
        else:
            # concave up
            RLcent = -RLfit[1] / (2*RLfit[0])

        # Dont jump fartther than max(sarr)
        if RLcent > np.max(sarr):
            RLcent = np.sign(RLcent) * np.max(sarr)


        if UBfit[0] < 0.0:
            # concave down
            print "UB concave down"
            UBcent = -2 * np.max(sarr) * np.sign(UBfit[1])
        else:
            # concave up
            UBcent = -UBfit[1] / (2*UBfit[0])

        # Dont jump fartther than max(sarr)
        if RLcent > np.max(sarr):
            UBcent = np.sign(UBcent) * np.max(sarr)

        # Rough conversion to pixels
        yguess += UBcent * binsize
        xguess -= RLcent * binsize

        print "{2:d}, update x0: {0}, y0: {1}".format(xguess, yguess, ii)
        if plotit:
            fig, ax = plt.subplots()
            ax.plot(sarr, UB, 'r', lw=1, label='UD')
            ax.plot(sarr, RL, 'b', lw=1, label='RL')
            ax.legend()
            plt.show()

        if np.sqrt((UBcent*binsize)**2 + (RLcent*binsize)**2) / binsize < 0.1:
            break

    return xguess, yguess


if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    # shot_number = 15676
    # fname = join(folder, "{0:07d}_000.nef".format(shot_number))
    # bg_fname = join(folder, "{0:07d}_001.nef".format(shot_number))
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    data = get_image_data(fname, bg_fname, color='b')

    ny, nx = data.shape
    x0 = nx/2.
    y0 = ny/2.
    x0, y0 = locate_center(data, x0, y0, plotit=True)
    binarr, ULsigarr, URsigarr, BLsigarr, BRsigarr = quick_ringsum(data, x0, y0, binsize=binsize)
    fig, ax = plt.subplots()
    ax.plot(binarr, ULsigarr, label='UL', lw=1)
    ax.plot(binarr, URsigarr, label='UR', lw=1)
    ax.plot(binarr, BLsigarr, label='BL', lw=1)
    ax.plot(binarr, BRsigarr, label='BR', lw=1)
    plt.legend(loc='upper left')
    plt.show(block=False)

    fig, ax = plt.subplots()
    i = ax.imshow(data, cmap='gray', origin="lower")
    fig.colorbar(i)
    ax.plot([0, data.shape[1]], [y0] * 2, 'b')
    ax.plot([x0] * 2, [0, data.shape[0]], 'r')
    ax.set_aspect(1.0)
    plt.show()

