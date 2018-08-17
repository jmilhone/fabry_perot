from __future__ import division, print_function
import concurrent.futures
from pymultinest import Analyzer
import numpy as np
from ..core import models
from ..tools import plotting, file_io
import matplotlib.pyplot as plt
import os.path as path


w0 = 468.619458
mu = 232.03806


def check_full_solver(finesse_folder):

    data_filename = path.join(finesse_folder, 'finesse_input.h5')

    data = file_io.h5_2_dict(data_filename)

    ix0 = data['fit_ix']['0']
    ix1 = data['fit_ix']['1']

    r0 = data['r'][ix0]
    sig0 = data['sig'][ix0]
    sig0_sd = data['sig_sd'][ix0]

    r1 = data['r'][ix1]
    sig1 = data['sig'][ix1]
    sig1_sd = data['sig_sd'][ix1]

    n_params = 6
    analyzer = Analyzer(n_params, outputfiles_basename=path.join(finesse_folder, "full_"))
    post = analyzer.get_equal_weighted_posterior()

    post = post[:, 0:-1]  # Don't include the last column, I forget what it means

    npost = post.shape[0]
    nsamples = 100
    samples = np.random.choice(npost, size=nsamples)

    sampled_post = post[samples, :]
    output0 = np.zeros((nsamples, len(r0)))
    output1 = np.zeros((nsamples, len(r1)))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures_map = dict()
        for idx, cube in enumerate(sampled_post):
            fut = executor.submit(model_wrapper, r0, r1, cube)
            futures_map[fut] = idx

        for future in concurrent.futures.as_completed(futures_map):
            idx = futures_map[future]

            val0, val1 = future.result()
            output0[idx, :] = val0
            output1[idx, :] = val1

    # calculate simple statistics for the fit
    output0_mean = np.mean(output0, axis=0)
    output0_low = np.min(output0, axis=0)
    output0_high = np.max(output0, axis=0)

    output1_mean = np.mean(output1, axis=0)
    output1_low = np.min(output1, axis=0)
    output1_high = np.amax(output1, axis=0)

    # Plotting fit results
    fig, ax = plt.subplots(figsize=(12,9))
    ax.fill_between(r0, output0_low, output0_high, color='C0', alpha=0.5)
    ax.plot(r0, output0_mean, color='C0', label='fit')

    ax.fill_between(r1, output1_low, output1_high, color='C0', alpha=0.5)
    ax.plot(r1, output1_mean, color='C0')

    ax.errorbar(r0, sig0, yerr=sig0_sd, color='C1', label='data')
    ax.errorbar(r1, sig1, yerr=sig1_sd, color='C1')
    ax.legend(fontsize=16)
    ax.set_xlabel('R (px)', fontsize=16)
    ax.set_ylabel('Counts', fontsize=16)
    ax.tick_params(labelsize=16)
    fig.tight_layout()
    plt.show()

    labels = ["L", "d", "F", "A", "B", "Ti"]
    fac = [0.004, 1.0, 1.0, 1.0, 1.0, 1.0]
    bins = ['auto', 'auto', 'auto', 'auto', 'auto', 'auto']
    for idx, label in enumerate(labels):
        print(label, ": ", np.nanmean(post[:, idx]*fac[idx]))
        fig, ax = plt.subplots()
        plotting.my_hist(ax, post[:, idx]*fac[idx], bins=bins[idx])
        ax.set_xlabel(label)
        plt.show(block=False)
    plt.show(block=True)


def model_wrapper(r0, r1, cube):
    vals0 = models.offset_forward_model(r0, cube[0], cube[1], cube[2], w0,
                                        mu, cube[3], cube[5], 0.0, coeff=0.05)
    vals1 = models.offset_forward_model(r1, cube[0], cube[1], cube[2], w0,
                                        mu, cube[4], cube[5], 0.0, coeff=0.05)
    return vals0, vals1
