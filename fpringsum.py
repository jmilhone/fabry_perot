import numpy as np
import matplotlib.pyplot as plt
import time
import rawpy
import scipy.optimize


def synthetic_data(plotit=0):
    x = np.arange(1, 1025, 1)  # synthetic data is a 1M CCD (like Andor camera).
    y = np.arange(1, 1025, 1)

    x0 = 500  # slightly off-center center
    y0 = 520

    X, Y = np.meshgrid(x - x0, y - y0)

    F = 30.  # finesse
    n = 1.  # index of refraction
    d = 0.88  # etalon spacing in mm
    lam = 468.6e-6  # wavelength of interest in mm (468.6e-6 for He and 488.0e-6 for Ar)
    L = 150.  # distance from lens to CCD (focal length of camera lens) in mm

    R = np.sqrt((X) ** 2. + (Y) ** 2.) * (13.3 / 1024.)  # conversion to mm from pixels

    A = (1. + ((4. * F ** 2.) / np.pi) * np.sin(
        2. * np.pi * ((n * d) / lam) * (L / np.sqrt(L ** 2. + R ** 2.))) ** 2.) ** -1.  # Ari function of rings

    if plotit:
        xmin = X[0, 0]
        xmax = X[0, -1]
        ymin = Y[0, 0]
        ymax = Y[-1, 0]

        plt.imshow(A, vmin=0, vmax=1, cmap='gray', extent=[xmin, xmax, ymin, ymax], interpolation='nearest',
                   origin='lower')
        plt.axis([xmin, xmax, ymin, ymax])
        plt.colorbar()
        plt.axes().set_aspect('equal')
        plt.show()

    return ((x0, y0), A)


def real_data(plotit=0):
    # rawpy is fantastic. It reads .nef's pretty fast. output_bps is the bit depth of the output (go with 16).
    # no_auto_bright=True turns off the auto-scaling of the image that makes the maximum pixel=2^16-1. This needs to
    # be set to True when doing background subtraction. adjust_maximum_thr=0. turns off some other weird normalization thing.

    rgb = rawpy.imread('Images/0014700_000.nef').postprocess(output_bps=16, no_auto_bright=True, adjust_maximum_thr=0.)

    center = (3066.33, 2031.42)  # center for 0014700_000.nef is (3066.33,2031.42) according to Cooper

    A = rgb[:, :, 2].astype(
        'float64')  # take only the blue pixel values, set as a float to remove worries about math later

    if plotit:
        plt.imshow(A, cmap='gray', interpolation='nearest', origin='lower')
        plt.colorbar()
        plt.axes().set_aspect('equal')
        plt.show()

    return (center, A)


@profile
def small_angle_sum(data, center, binsize=1, plotit=0, printit=0):
    data = data / data.max()  # Normalize the data to one

    t1 = time.time()  # Timing function is just to show off

    nx = data.shape[1]  # get the size and shape of 'data'
    ny = data.shape[0]
    npts = data.size

    x = np.arange(1, nx + 1, 1)  # x and y arrays of the pixel locations on CCD
    y = np.arange(1, ny + 1, 1)

    x0, y0 = center  # center must be provided in this version

    X, Y = np.meshgrid(x - x0, y - y0)  # X and Y are arrays same shape and size as data, but listing the
    # corresponding x,y pixel location values wrt center

    del x, y  # these are fairly large arrays that we don't need anymore so this removes them from memory

    xmin = X[0, 0]
    xmax = X[0, -1]  # xmin,xmax,ymin,ymax are the extents of the chip wrt center
    ymin = Y[0, 0]
    ymax = Y[-1, 0]

    R = np.sqrt(
        X ** 2. + Y ** 2.)  # R is the same shape and size as data, but lists the R location of pixels wrt center

    del X, Y  # these were big so deleting them from memory is nice

    # plt.imshow(R,vmin=0,vmax=R.max(),extent=[xmin,xmax,ymin,ymax],interpolation='nearest',origin='lower')
    # plt.colorbar()
    # plt.axes().set_aspect('equal')
    # plt.show()

    rmax = np.min(
        np.abs([xmin, xmax, ymin, ymax]))  # the maximum r we want to go to is the nearest edge of CCD to center

    # The small angle approximation calls for bins in radius of equal area. Thus binsize of r decreases as a function of r.
    # In order to maintain equal area we require that binsize_r(r)~sqrt(c+r^2)-r where c is a constant. The smallest binsize
    # we will allow is the binsize keyword supplied (default to 1) and this binsize will occur at r=rmax. This tells us that
    # c=2*rmax*binsize+binsize^2 and, therefore, binsize_r(r)=sqrt(r^2+2*rmax*binsize+binsize^2)-r. Now we can turn this into an
    # iterable function such that binarr(i)=sqrt(2*i*rmax*binsize+i*binsize^2) with i->0-imax. Imax can be found from solving for
    # i with binarr=rmax. In python, indices start at 0 so the formula used in the below function has i->i+1.

    imax = int((rmax ** 2 - 2 * rmax - 1) / (1 + 2 * rmax) / binsize)
    binarr = np.fromfunction(lambda i: np.sqrt(2. * (i + 1.) * rmax * binsize + (i + 1.) * binsize ** 2.), (imax,),
                             dtype='float64')

    # binarr is the array of equal area bins in r-space that is proportional to lambda space in the small angle approximation

    R = R.reshape((npts,))  # ok now we are going to flatten both R and data into 1d arrays
    data = data.reshape(
        (npts,))  # since the flatten is the same method and they were the same shape, they should still match
    rarr_ix = np.argsort(R,
                         kind='mergesort')  # this is a pretty slow line but we need to sort R from smallest to largest
    R = R[rarr_ix]
    data = data[rarr_ix]

    R = R[R < rmax]  # trim off anything greater than the rmax of the bins we made
    data = data[0:R.size]  # do the same correspondingly to data

    # This is the actual ring summing (in one line, I might add). We add a 0 to the beginning of binarr so we have a bottom
    # bin. Histogram bins R into bins given by [0,binarr] but instead of counting the number of items from R in each bin, it
    # assigns each of them a weight given by data and sums all the bins. Exactly like we want for a ring sum.

    sigarr, _ = np.histogram(R, bins=np.concatenate((np.array([0.]), binarr)), weights=data)

    t2 = time.time()  # just to keep track of how fast this is
    if printit:
        print 'Summed a ring in ' + np.str(np.round((t2 - t1) / 0.001) * 0.001) + ' s'

    if plotit:
        plt.plot(binarr, sigarr, 'b')
        plt.show()

    return (binarr, sigarr)


def center_minimize(center, data):
    x0 = int(center[0])
    y0 = int(center[1])

    ym = np.sum(np.abs(data[y0::, :] - data[:-y0, :][::-1, :]))

    xm = np.sum(np.abs(data[:, x0::] - data[::, :-x0][:, ::-1]))

    return xm + ym


def find_center(data):
    (ny, nx) = data.shape

    (x0, y0) = (int(nx / 2), int(ny / 2))  # start with center of chip as a guess for the center


(center, data) = real_data()
# (center,data)=synthetic_data()
(binarr, sigarr) = small_angle_sum(data, center, binsize=0.1)
