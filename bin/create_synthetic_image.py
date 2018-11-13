from __future__ import division, print_function
#import sys
#sys.path.append("../")
from fabry.plasma import plasma
from fabry.core import models, ringsum
from fabry.tools import file_io, images
import numpy as np
import concurrent.futures
import matplotlib.pyplot as plt
import time
from numba import jit

def main(w, mass, temp, v_outer, r_outer, impact_factor, Lnu, fname, rmax=42.0):
    warr, spectrum = plasma.calculate_pcx_chord_emission(impact_factor, temp, w, mass, Lnu, v_outer, rmax=rmax, nr=500,
                                                         nlambda=2000, Lne=4.0, R_outer=r_outer)

    L = 150.0 / 0.004
    d = 0.88
    F = 20.7
    
    print('Impact Factor', impact_factor)
    print('Ti: ', temp)
    print('w0: ', w)
    print('mu: ', mass)
    print('Lnu: ', Lnu)
    print('Vouter: ', v_outer)
    print('Router: ', r_outer)
    print('Rmax: ', rmax)
    print('L: ', L)
    print('d: ', d)
    print('F: ', F)
    
    r, theta, x = plasma.calculate_r_theta_x_from_impact_factor(impact_factor, rmax=rmax, npts=500)
    vel = plasma.pcx_velocity_profile(r, Lnu, r_outer, v_outer)
    plt.plot(r, vel)
    plt.show()

    raw_input('pausing...')
    nx = 6000
    ny = 4000

    data = np.zeros((ny, nx))
    x = np.arange(1, nx+1, 1)
    y = np.arange(1, ny+1, 1)

    x0 = 3000
    y0 = 2000

    XX, YY = np.meshgrid(x, y)
    R = np.sqrt((XX-x0)**2 + (YY-y0)**2)

    nprocs = 32 

    split_r = np.array_split(R.flatten(), nprocs)
    results = [0.0 for _ in split_r]
    t0 = time.time()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        fut_to_spl = dict()
        for idx, split in enumerate(split_r):
             future = executor.submit(models.general_model, split, L, d, F, warr, spectrum)
             fut_to_spl[future] = idx
        for future in concurrent.futures.as_completed(fut_to_spl):
             current_idx = fut_to_spl[future]
             output = future.result()
             results[current_idx] = output

    results = np.concatenate(results)
    results.shape = (ny, nx)

    results *= 1000.0 / results.max()

    t1 = time.time()
    print(t1 - t0, 'seconds')
    fig, ax = plt.subplots()
    im = ax.imshow(results)
    plt.colorbar(im, ax=ax)
    plt.show()

    x0, y0 = ringsum.locate_center(results, xguess=3000, yguess=2000)
    rr, ss, _ = ringsum.ringsum(results, x0, y0)

    rrr = np.linspace(1, 1750, 5000)
    actual = models.general_model(rrr, L, d, F, warr, spectrum)

    fig, ax = plt.subplots()
    ax.plot(rr, ss/ss.max(), 'C0', label='ringsum')
    ax.plot(rrr, actual/actual.max(), 'C1', label='actual')
    ax.legend()
    plt.show()

    #results = results.astype(int)
    #results = np.asarray(results, dtype=np.int16)
    image_data = {'image': results, 'L': L, 'd': d, 'F': F, 'Vouter': v_outer, 
            'R_outer': r_outer, 'Ti': temp, 'w0': w, 'mu': mass, 'Lnu': Lnu, 
            'impact_factor': impact_factor
            }
    file_io.dict_2_h5(fname, image_data) 


def add_noise(input_data):
    output_data = np.zeros_like(input_data)
    n, m = input_data.shape
    noise_scale = np.sqrt(input_data)
    for i in range(n):
        print(i)
        for j in range(m):
            output_data[i, j] = input_data[i, j] + np.random.normal(loc=0.0, scale=noise_scale[i, j], size=1)

    return output_data # output_data.astype(np.int16)

if __name__ == "__main__":
    w0 = 487.98634
    mu = 39.948
    Ti = 0.5
    Vouter = 1000.0
    Router = 35.0
    Rmax = 42.0
    impact_fac = 35.0
    ne = 2e17
    nn = 8e17
    # mom_dif_length = 20.0
    mom_dif_length = 100.0*plasma.Lnu(ne*nn, Ti, mu=40.0)
    fname = 'pcx_vel_profile_{0:2.0f}.h5'.format(impact_fac)
    print(fname)

    main(w0, mu, Ti, Vouter, Router, impact_fac, mom_dif_length, fname, rmax=Rmax)
    data = file_io.h5_2_dict(fname)
    data['image'] = add_noise(data['image'])

    file_io.dict_2_h5(fname, data)

    # fig, ax = plt.subplots()
    # im = ax.imshow(data['image'])
    # plt.colorbar(im)
    # plt.show()

