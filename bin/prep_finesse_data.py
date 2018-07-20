from __future__ import print_function, division, absolute_import
import numpy as np
import argparse
from fabry.tools.file_io import h5_2_dict, dict_2_h5
from fabry.tools.plotting import ringsum_click
import matplotlib.pyplot as plt
from os.path import join, abspath


def get_finesse_region(r_array, sig_array, sig_err_array, plot_fit_region=True):
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
    index1 = np.abs(edge_r[0] - r_array).argmin()
    index2 = np.abs(edge_r[1] - r_array).argmin()
    fit_indices = np.arange(index1, index2+1, 1, dtype=int)

    if plot_fit_region:
        fig, ax = plt.subplots()
        ax.errorbar(r_array, sig_array, yerr=sig_err_array , color='C0', label='Signal', ecolor='C2', errorevery=5)
        ax.plot(r_array[fit_indices], sig_array[fit_indices], 'C1', label='Signal to Fit', zorder=100)
        ax.axhline(sig_array.max() * 0.15 / (1.0 + 21.0), color='k')
        ax.legend(frameon=False)
        plt.show()
    return fit_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Determine finesse fit region and write to hdf5 file')
    parser.add_argument("folder", type=str, help='Folder containing output (ringsum.h5) of process_image.py')
    args = parser.parse_args()

    folder = abspath(args.folder)
    fname = join(folder, 'ringsum.h5')

    data = h5_2_dict(fname)
    print(data.keys())
    r = data['r']
    sig = data['sig']
    data['sig_sd'] = np.sqrt(data['sig_sd']**2 + (0.01*sig)**2) # this is the error contribution from the center error
    min_loc = np.argmin(sig)
    # data['sig'] -= sig.min()

    # data['sig_sd'] = np.sqrt(data['sig_sd']**2 + data['sig_sd'][min_loc]**2) # Adding in the offset error

    ix = get_finesse_region(r, sig, data['sig_sd'], plot_fit_region=True)

    data['fit_ix'] = ix

    fig, ax = plt.subplots()
    ax.plot(r, data['sig_sd'] / data['sig'])
    plt.show()

    dict_2_h5(join(folder, 'finesse_input.h5'), data)

