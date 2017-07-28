from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
import argparse
#from He_model import main as he_mod
from scipy.optimize import differential_evolution
import He_model2 as model

#L = 149.496800727 / .004
#d = 0.88342569509 

#L = 37379.4286
#d = .883669660669

#L = 149.673848949 / .004 
#d = 0.883181102804 

#L = 149.688242473  / .004
#d = 0.883425147574 

#L = 149.69264395 / .004 
#d = 0.883425109538 

#L = 149.439432006  / .004
#d = 0.883426095104

#L = 148.893784919 / .004
#d = 0.878060370172


#L = 148.642754807/.004 
#d = 0.872202304089


def chi_sq(a, xx, yy):
    #L = 149.496800727 / .004
    #d = 0.88342569509 

    #L = 148.717940474 / .004
    #d = 0.872201537841
    L = 148.715306908/ .004
    d = 0.872201921125
    #L = 148.599558491 / .004
    #d = 0.872203104696
    doffset = a[0]
    d += doffset
    Ti = a[1]
    amplitudes = a[2:]

    ymodel = model.model_output(xx, L, d, Ti, amplitudes)
    
    return np.sum((yy - ymodel)**2 / 0.05**2)


def calculate_peaks(L, d, w, norders=4):
    m = 2.e6 * d / w
    m0 = np.floor(m)
    return [L * np.sqrt( m**2 / (m0 - j)**2 - 1.0) for j in range(norders)]


def main(fname, L, d, bgname=None, x0=None, y0=None, binsize=1.0, subtract_min=False):
    data = im.get_image_data(fname, bgname=None, color='b')
    
    if x0 is None:
        x0 = data.shape[1]/2
    if y0 is None:
        y0 = data.shape[0]/2

    #im.quick_plot(data)
    #plt.show()
    #x0 = 3043.988356 
    #y0 = 2001.86956396
    #x0 = 3036.27977126 
    #y0 = 2023.21834704
    #x0 = 3042.23040935
    #y0 = 2003.94310613
    x0 = 3043.95786886
    y0 = 2001.8908683
    #x0, y0 = rs.locate_center(data, x0, y0, binsize=0.1, plotit=False)
    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)

    if bgname:
        bgdata = im.get_image_data(bgname, bgname=None, color='b')
        _, bg_sig = rs.quick_ringsum(bgdata, x0, y0, binsize=binsize, quadrants=False)
        sigarr -= bg_sig

    #im.quick_plot(data-bgdata)
    #ax = plt.gca()
    #ax.axhline(y0)
    #ax.axvline(x0)
    #plt.show()
    if subtract_min:
        sigarr -= sigarr.min()

    #V = 5.0
    V = 0.0
    print V
    Ti = 0.5

    amplitudes = [0.32, 0.2, 1.13, 1.61, 1.0]

    #rmin = 670
    rmin = 650
    rmax = 970
    imin = np.abs(binarr - 670).argmin()
    imax = np.abs(binarr - 970).argmin()

    # Normalize data for fitting
    x = binarr[imin:imax+1].copy()
    y = sigarr[imin:imax+1].copy()
    ymax = y.max()
    y /= ymax

    bounds = [(-5e-5, 7.5e-5), (0.025, 3.0), (0.0,3.0), (0.0,3.0), (0.0,3.0), 
            #(0.0,3.0) ] 
            (0.0,3.0), (0.0,3.0)] 

    res = differential_evolution(chi_sq, bounds, args=(x, y))
    print res

    amplitudes = res['x'][2:]
    #d = d + res['x'][0]
    Ti = res['x'][1]
    print len(res['x'])
    print "Ti: ",Ti
    print "d: ", d
    print amplitudes

    print L*.004, d
    ymodel = model.model_output(x, L, d, Ti, amplitudes, F=20.7)
    #plt.plot(x, y)
    plt.plot(binarr, sigarr, 'g')
    plt.plot(x, ymax*ymodel, 'r')
    plt.show()

    #rr, ss = he_mod(L, d+4.3e-5, V=V, Ti=Ti)
    #rr = np.linspace(0, 2000, 10000)
    #amplitudes = [0.32, 0.2, 1.13, 1.61, 1.0]
    #ss = model.model_output(rr, L, d+4.3e-5, Ti, amplitudes)
    #ss *= np.exp(-(rr / 2400)**2 )
    #fig, ax = plt.subplots()
    #ax.plot(binarr, sigarr, 'b')
    #ax1 = ax.twinx()
    #ax1.plot(rr, ss, 'r')
    #ax.plot(binarr[imin:imax], sigarr[imin:imax], 'g')
    #ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    #lo, hi = ax.get_ylim()
    #ax.set_ylim(0, hi)
    #plt.show()


if __name__ == "__main__":
    w0 = [468.564736669, 465.4559]
    #w0 = [468.6195, 467.3660927]
    #for w in w0:
    #    print "\n", w, "nm"
    #    pks = calculate_peaks(L, d, w, norders=4)
    #    print [x**2 for x in pks]
    #    #print pks
    fname = "0009174_000.nef"
    bgname = "0009174_001.nef"
    #fname = "0009334_000.nef"
    #bgname = "0009334_001.nef"
    #fname = "0009370_000.nef"
    #bgname = "0009370_001.nef"
    #fname = "0009395_000.nef"
    #bgname = "0009395_001.nef"

    #fname = "0010641_000.nef"
    #bgname = "0010641_001.nef"
    #fname = "Th_lamp_468_calib_5m.nef"
    
    folder = "Images/"
    fname = join(folder, fname)
    bgname = join(folder, bgname)
    #bgname = None
    xguess = 3043.0
    yguess = 2001.0
    #L = 148.717940474 / .004
    #d = 0.872201537841
    L = 148.715306908/ .004
    d = 0.872201921125
    #L = 148.599558491 / .004
    #d = 0.872203104696

    #L = 149.496800727 / .004
    #d = 0.88342569509 
    main(fname, L, d, bgname=bgname, x0=xguess, y0=yguess, binsize=2.0, subtract_min=True)




