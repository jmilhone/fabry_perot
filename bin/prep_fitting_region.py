from __future__ import print_function, division, absolute_import
#import sys
#sys.path.append("../")
import numpy as np
import argparse
from fabry.tools.file_io import h5_2_dict, dict_2_h5
from fabry.tools.plotting import ringsum_click
import matplotlib.pyplot as plt
from os.path import join, abspath
from scipy.special import wofz
from scipy.optimize import curve_fit


def voigt(x, x0, sigma, gamma):
    z = (x-x0 + 1j*gamma) / sigma / np.sqrt(2) 
    w = wofz(z)
    return np.real(w) / sigma / np.sqrt(2*np.pi)

def double_voigt(x, x1, x2, sigma1, sigma2, gamma, amp1, amp2):
    return voigt(x, x1, sigma1, gamma)*amp1 + voigt(x,x2,sigma2,gamma)*amp2

def remove_voigt(x, y):
    edge_r, s = ringsum_click(x, y, title='Click to the two peaks to fit')

    a0 = [edge_r[0], edge_r[1], 1.0, 1.0, 1.0, s[0], s[1]]
    bounds = ([edge_r[0]-10, edge_r[1]-20, 0.01, 0.01, 0.001, 0.5*s[0], 0.1*s[1]], 
              [edge_r[1]+10, edge_r[1]+20, 1e2, 1e2, 1e4, 10*s[0], 10*s[1]])
    print(bounds)
    popt, pcov = curve_fit(double_voigt, x, y, p0=a0, 
            bounds=bounds)
    print(popt)

    fig, ax = plt.subplots()
    ax.plot(x, y, 'b')
    ax.plot(x, double_voigt(x, *popt),'r')
    plt.show()



def get_fitting_region(r_array, sig_array, sig_err_array, plot_fit_region=True):
    """
    User clicks the left and right edges of the region they would like to fit for the finesse solver

    Arguments:
        r_array (np.ndarray): radius array
        s_array (np.ndarray): signal array
        plot_fit_region (bool): True if the user wants to see the fitted region overplotted the signal.

    Returns:
        (np.ndarray): Array of indices for the fitting region
    """
    edge_r, _ = ringsum_click(r_array, sig_array, title='Click to the left and right of the finesse fitting region')
    indices = [np.abs(edge-r_array).argmin() for edge in edge_r]

    left = indices[0::2]
    right = indices[1::2]

    fit_indices = [range(L, R+1) for (L, R) in zip(left, right)]
    # Quick modification
    fit_indices = np.concatenate(fit_indices)
    fit_indices_dict = {'0': fit_indices}
    #fit_indices_dict = {str(idx): range(L, R+1) for (idx, (L, R)) in enumerate(zip(left, right))}

    #index1 = np.abs(edge_r[0] - r_array).argmin()
    #index2 = np.abs(edge_r[1] - r_array).argmin()
    #fit_indices = np.arange(index1, index2+1, 1, dtype=int)

    if False:#plot_fit_region:
        fig, ax = plt.subplots()
        ax.errorbar(r_array, sig_array, yerr=sig_err_array, color='C0', label='Signal', ecolor='C2', errorevery=5)
        #ax.plot(r_array[fit_indices], sig_array[fit_indices], 'C1', label='Signal to Fit', zorder=100)
        for slice in fit_indices:
            ax.axvspan(r_array[slice][0], r_array[slice][-1], color='C2', alpha=0.5)
        #ax.axhline(sig_array.max() * 0.5 / (1.0 + (2*21.0 / np.pi)**2), color='k')
        ax.legend(frameon=False)
        plt.show()
    return fit_indices_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Determine finesse fit region and write to hdf5 file')
    parser.add_argument("folder", type=str, help='Folder containing output (ringsum.h5) of process_image.py')
    parser.add_argument("filename", type=str, help="Filename to be written in <folder>/")
    args = parser.parse_args()

    folder = abspath(args.folder)
    fname = join(folder, 'ringsum.h5')

    data = h5_2_dict(fname)
    print(data.keys())
    r = data['r']
    sig = data['sig']
    #data['sig_sd'] = np.sqrt(data['sig_sd']**2 + (0.02*sig)**2) # this is the error contribution from the center error
    data['sig_sd'] = np.sqrt(data['sig_sd']**2 + (0.005*sig)**2) # this is the error contribution from the center error
    
    #ix = get_fitting_region(r, sig, data['sig_sd'], plot_fit_region=True)['0']
    #print(ix)
    #remove_voigt(r[ix], sig[ix])

    #print("***********************************************") 
    #print("Adding error to the region between ThI and ArII")
    #print("***********************************************") 
    #r1 = 540
    #r2 = 600
    #i1 = np.abs(r - r1).argmin()
    #i2 = np.abs(r - r2).argmin()+1
    #sl = slice(i1, i2)
    #data['sig_sd'][sl] = np.sqrt(data['sig_sd'][sl]**2 - (0.01*data['sig'][sl])**2 + (0.02*data['sig'][sl])**2)

    #min_loc = np.argmin(sig)
    # data['sig'] -= sig.min()
    # data['sig'] -= 0.8 * sig.min()
    # data['sig_sd'] = np.sqrt(data['sig_sd']**2 + data['sig_sd'][min_loc]**2) # Adding in the offset error

    ix = get_fitting_region(r, sig, data['sig_sd'], plot_fit_region=True)

    data['fit_ix'] = ix

    fig, ax = plt.subplots()
    ax.plot(r**2, data['sig'])
    ax.plot(r[ix['0']]**2, data['sig'][ix['0']])
    ax.fill_between(r**2, data['sig'] - data['sig_sd'], data['sig'] + data['sig_sd'], alpha=0.5)
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(r[ix['0']], data['sig_sd'][ix['0']] / data['sig'][ix['0']])
    plt.show()

    #dict_2_h5(join(folder, 'finesse_input.h5'), data)
    # dict_2_h5(join(folder, 'test_finesse_input.h5'), data)
    dict_2_h5(join(folder, args.filename), data)

