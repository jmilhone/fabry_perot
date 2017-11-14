from __future__ import division
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import argparse
import image_helpers as im
import ring_sum as rs
import fitting
from os.path import join
import time
import plottingtools.core as ptools
import json
from scipy.stats import norm
import fp_helpers

class ListBuilder():
    def __init__(self, fig):
        self.x = []
        self.fig = fig
        self.cid = fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        self.x.append(event.xdata)


def determine_fit_range(x, y, L, R, thres=0.15):
    iL = np.abs(x - L).argmin()
    iR = np.abs(x - R).argmin()
    indices = range(iL,iR+1)
    i0 = np.argmax(y[indices]) + indices[0]
    #plt.plot(x, y)
    #plt.plot(x[indices], y[indices])
    #plt.plot(x[i0], y[i0], 'or')
    #plt.axhline(thres*y[i0], color='k')
    #plt.show()
    LL = np.max(np.where(y[0:i0] < thres*y[i0]))
    RR = np.min(np.where(y[i0:] < thres*y[i0])) + i0
    return LL, RR

def find_maximum(x, y):
    # I am assuming there is only going to be one minimum!
    dy = np.diff(y)
    xx = 0.5*(x[0:-1] + x[1:])
    i = np.where(np.logical_and(dy[0:-1]>0.0, dy[1:]<0.0))[0]
    if len(i) > 1:
        #print "adjust i to be one value"
        vals = y[i]
        idx = np.argmax(vals)
        i = i[idx]
    slope = (dy[i+1] - dy[i]) / (xx[i+1] - xx[i])

    if slope < 0.0:
        pk = xx[i] - dy[i] / slope
        #return pk[0]
        return pk
    else:
        return None



def argon_calib(fname, x0, y0, bg=None, plotit=True, color='b', binsize=0.5, write_peaks=False, folder=None):
    if fname[-3:].lower() == 'nef':
        data = im.get_image_data(fname, color=color)
    else:
        data = np.load(fname)

    if plotit:
        fig, ax = plt.subplots()
        ax.set_ylim(0, 4000)
        ax.set_xlim(0, 6000)
        ax.axis('equal')
        cb = ax.imshow(data, cmap='gray')
        plt.colorbar(cb, ax=ax)
        plt.show()

    x0, y0 = rs.locate_center(data, x0, y0, plotit=False)

    ny, nx = data.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)

    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)
    theta = np.arctan2(yy, xx)

    theta_min = np.deg2rad(120.0) #np.arctan2(-2730.0+y0, 2364.0-x0)
    theta_max = np.deg2rad(150.0) #np.arctan2(-2987.0+y0, 2925.0-x0)
    
    ii = np.where(np.logical_and(theta_min < theta, theta_max > theta))
    # data[ii] = 0.0
    if plotit:
        fig, ax = plt.subplots()
        ax.set_ylim(0, 4000)
        ax.set_xlim(0, 6000)
        ax.axis('equal')

        cb = ax.imshow(data, cmap='gray')
        ax.axhline(y0, color='b')
        ax.axvline(x0, color='r')
        plt.colorbar(cb, ax=ax)
        plt.show()

    # compute ring sum and use center of bins for r arrray
    binarr, sig = rs.quick_ringsum(data, x0, y0, binsize=1.0, quadrants=False)
    r = np.concatenate(([0.0], binarr))
    rdiff = np.diff(r)
    r = 0.5 * (r[0:-1] + r[1:])

    # Handle background image here
    if bg:
        if bg[-3:].lower() == 'nef':
            bgdata = im.get_image_data(bg, color=color)
        else:
            bgdata = np.load(bg)
        _, bg_sig = rs.quick_ringsum(bgdata, x0, y0, binsize=1.0, quadrants=False)
        sig -= bg_sig

    boundaries = {'Ar left': [],
                  'Ar right': [],
                  'Th left': [],
                  'Th right': []}
    fig_names = {'Ar left': "Click to the left of the Ar peaks",
                 'Ar right': "Click to the right of the Ar peaks",
                 'Th left': "Click to the left of the Th peaks",
                 'Th right': "Click to the right of the Th peaks"}

    fig, ax = plt.subplots()
    ax.plot(r, sig, 'b')
    #ax.plot(r**2, sig, 'b')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_xlabel("R (px)")
    #ax.set_xlabel(r"$R^2$ (px${}^2$)")
    ax.set_ylabel("Counts")
    plt.show()

    fname = "calibration_peaks_location_data_2017_10_30.json"
    w0 = ['487.873302', '487.98634']
    w = [float(x) for x in w0]

    if write_peaks:
        for x in boundaries:
            fig, ax = plt.subplots()    
            ax.plot(r, sig, 'b')
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.set_xlabel("R (px)")
            ax.set_ylabel("Counts")
            boundary_builder = ListBuilder(fig)
            plt.title(fig_names[x])
            plt.show()

            fig.canvas.mpl_disconnect(boundary_builder.cid)
            boundaries[x] = boundary_builder.x
        print boundaries
        peak_location_data = {"left":{}, "right":{}}
        peak_location_data['left']['487.873302'] = boundaries['Th left'] 
        peak_location_data['right']['487.873302'] = boundaries['Th right'] 
        peak_location_data['left']['487.98634'] = boundaries['Ar left'] 
        peak_location_data['right']['487.98634'] = boundaries['Ar right'] 
        with open(fname, 'w') as outfile:
            json.dump(peak_location_data, outfile, indent=4, separators=(',', ': '))
    else:
        with open(fname, 'r') as infile:
            peak_location_data = json.load(infile)


    norders = len(peak_location_data['left'][w0[0]])
    ar_peaks = []
    th_peaks = []
    for i in xrange(norders):
        L_Ar, R_Ar = determine_fit_range(r, sig, peak_location_data['left'][w0[1]][i],
                peak_location_data['right'][w0[1]][i], thres=0.4)
        L_Th, R_Th = determine_fit_range(r, sig, peak_location_data['left'][w0[0]][i],
               peak_location_data['right'][w0[0]][i], thres=0.6)

        #fig, ax = plt.subplots()
        #ax.plot(r, sig)
        #ax.axvline(peak_location_data['left'][w0[0]][i], color='m')
        #ax.axvline(peak_location_data['right'][w0[0]][i], color='c')
        #ax.axvline(r[L_Th], color='k')
        #ax.axvline(r[R_Th], color='k')
        #plt.show()

        # This is a stupid fix...
        temp = find_maximum(r[L_Th:R_Th], sig[L_Th:R_Th])
        if isinstance(temp, np.ndarray):
            temp = list(temp)
            th_peaks += temp
        else:
            th_peaks.append(temp)

        temp = find_maximum(r[L_Ar:R_Ar], sig[L_Ar:R_Ar])
        if isinstance(temp, np.ndarray):
            temp = list(temp)
            ar_peaks += temp
        else:
            ar_peaks.append(temp)

    print ar_peaks
    print th_peaks
    
    """
    Run a new ring sum with the specified bin size.  The small bin size leads 
    to noisy data but we want that for multinest, but not for finding minimum 
    of derivatives.
    """
    binarr, sig = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    r = np.concatenate(([0.0], binarr))
    rdiff = np.diff(r)
    r = 0.5 * (r[0:-1] + r[1:])

    # Handle background image here
    if bg:
        _, bg_sig = rs.quick_ringsum(bgdata, x0, y0, binsize=binsize, quadrants=False)
        sig -= bg_sig

    # print "Arbitrary offset subtraction"
    # sig -= 25000.0
    # sig -= 2.5e6

    fit_range = []
    for i in xrange(norders):
        L_Ar, R_Ar = determine_fit_range(r, sig, peak_location_data['left'][w0[1]][i],
                peak_location_data['right'][w0[1]][i], thres=0.3)
        L_Th, R_Th = determine_fit_range(r, sig, peak_location_data['left'][w0[0]][i],
               peak_location_data['right'][w0[0]][i], thres=0.3)
        fit_range.append(range(L_Th, R_Ar+1))



    fig, ax = plt.subplots()
    #ax.plot(r**2, sig)
    ax.plot(r, sig)
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    for i in xrange(norders):
        #ax.plot(r[fit_range[i]]**2, sig[fit_range[i]])
        #ax.axvline(ar_peaks[i]**2, color='c')
        #ax.axvline(th_peaks[i]**2, color='m')
        ax.plot(r[fit_range[i]], sig[fit_range[i]])
        ax.axvline(ar_peaks[i], color='c')
        ax.axvline(th_peaks[i], color='m')
    # ax.set_xlabel(r"R${}^2$ (px${}^2$)")
    ax.set_xlabel(r"R (px)")
    ax.set_ylabel("Counts")
    if folder:
        #ax.set_xlim(0.0, 3.e6)
        savefolder = join(folder, "Plots/")
        fp_helpers.make_directory(savefolder)
        fig.savefig(join(savefolder, "spectra_fit_range_peaks.pdf"))
    plt.show()

    return r, sig, ar_peaks, th_peaks, fit_range


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Reads calibration image and performs initial processing for MutliNest solvers",
            )
    parser.add_argument("--finesse-dir", "-f", action='store', type=str, dest='finesse_dir', 
            help="Name of finesse directory located in Calibration/Finesse_saves/")
    parser.add_argument("--Ld-dir", "-d", action='store', type=str, dest='Ld_dir', 
            help="Name of Ld directory located in Calibration/Ld_saves/")
    parser.add_argument("-x", action='store', type=float, help="Initial guess for the x center")
    parser.add_argument("-y", action='store', type=float, help="Initial guess for the y center")
    parser.add_argument("--argon-calib-fname", "-A", action='store', type=str, dest='argon_fname')
    parser.add_argument("--argon-calib-bg-fname", action='store', type=str, default=None, 
            dest="argon_bg_fname")
    args = parser.parse_args()

    finesse_dir = join("Calibration/Finesse_saves/", args.finesse_dir, "")
    Ld_dir = join("Calibration/Ld_saves/", args.Ld_dir, "")

    fp_helpers.make_directory(finesse_dir)
    fp_helpers.make_directory(Ld_dir)

    print args.argon_fname
    binsize = 0.1 
    rarr, counts, Ar_peaks, Th_peaks, fit_indices = argon_calib(
            args.argon_fname, args.x, args.y, bg=args.argon_bg_fname, plotit=True, binsize=binsize, folder=finesse_dir)

    peak_data = {'Ar': Ar_peaks,
                 'Th': Th_peaks,}

    calib_data = {'r': rarr.tolist(),
                  'sig': counts.tolist(),
                  'idx': fit_indices}

    with open(join(Ld_dir, "calibration_peaks.json"), 'w') as peakfile:
        json.dump(peak_data, peakfile, indent=4, separators=(',', ': ')) 

    with open(join(finesse_dir, "calib_data.json"), 'w') as calibfile:
        json.dump(calib_data, calibfile, indent=4)



