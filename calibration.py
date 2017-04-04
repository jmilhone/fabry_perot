import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import minimize_scalar
import time

def find_initial_Lguess(r1, r2, wavelength, d):
    f = lambda x: np.abs(2.*d*1.e6*x/wavelength * (1./np.sqrt(x**2 + r2**2) - 1./np.sqrt(x**2+r1**2)) + 1.0)
    # (d in mm)*1.e6 = d in nm

    opt = minimize_scalar(f, bounds=(1e4, 1e6), method='bounded')

    # L = np.logspace(4, 6, 10000)
    # plt.plot(L, f(L))
    # plt.xlim(10000, 50000)
    # plt.show()
    # print opt
    # print opt.x
    return opt.x

if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    data = im.get_image_data(fname, bg_fname, color='b')

    #Very close guess for center
    # x0, y0 = (3069., 2032.)
    # REALLY close
    x0, y0 = (3069.68678585, 2032.85668627)

    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=False)
    duc = 0.88
    c_lambda = 487.8733
    r1 = peaks[0]
    r2 = peaks[1]
    Luc = find_initial_Lguess(peaks[0], peaks[1], c_lambda, duc)
    print Luc
    m1 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc**2 + r1**2) / c_lambda
    m2 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc**2 + r2**2) / c_lambda
    print m1, m2

    m1 = round(m1)
    m2 = round(m2)
    assert m1 - m2 == 1, "The rounding of m1 and m2 in the initial calibration failed"

    # Recalculate L and d with m1 and m2 being integers
    L = np.sqrt(((m2*r2)**2 - (m1*r1)**2)/(m1**2 - m2**2))
    print L
    d = m1 * c_lambda * np.sqrt(L**2 + r1**2) / (2*L*1.e6)
    print d


    t1 = time.time()
    ny, nx = data.shape
    x = np.arange(1, nx+1, 1)
    y = np.arange(1, ny+1, 1)
    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)
    delta_lambda = 0.0001  # in nm
    n_dlambda = 512  # This is Cooper's choice
    c_lambda_arr = np.arange(-n_dlambda/2, n_dlambda/2, 1)*delta_lambda + c_lambda
    lambda_min = c_lambda_arr.min()
    lambda_max = c_lambda_arr.max()
    npeaks = len(peaks)
    R = R.flatten()
    data = data.flatten()
    print time.time() - t1
    for i in range(npeaks):
        t2 = time.time()
        m = m1 - i
        rmin = np.sqrt((2*d*L*1.e6/(m*lambda_max))**2 - L**2)
        rmax = np.sqrt((2*d*L*1.e6/(m*lambda_min))**2 - L**2)
        inds = np.where(np.logical_and(R > rmin, R<rmax))[0]

        r0 = 2.*d*1.e6*L / m / delta_lambda
        # print r0, rmin, rmax

        dr = np.sqrt((1./(1./np.sqrt(L**2 + rmin**2) - 1./r0))**2 - L**2) - rmin
        rR = [rmin, rmin + dr]
        j = 0
        while rR[-1] < rmax:
            j += 1
            dr = np.sqrt((1. / (1. / np.sqrt(L ** 2 + rR[-1] ** 2) - 1. / r0)) ** 2 - L ** 2) - rR[-1]
            # print j, dr, rR[-1], rmax
            rR += [rR[-1] + dr]

        # print len(rR), len(c_lambda_arr)

        sig, _ = np.histogram(R[inds], bins=rR, weights=data[inds])
        # to correspond to my calibration wavelength array, sig is actually backwards...
        sig = sig[::-1]
        print time.time() - t2
        plt.plot(c_lambda_arr, sig)
        plt.title("idx: {0} order: {1}".format(i, m))
        plt.show()

