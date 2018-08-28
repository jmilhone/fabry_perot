from __future__ import division, print_function
import sys
sys.path.append("../")
from fabry.plasma import plasma
from fabry.core import models
import numpy as np
import concurrent.futures
import matplotlib.pyplot as plt

def main(w, mass, temp, v_outer, r_outer, impact_factor, Lnu, rmax=42.0):
    warr, spectrum = plasma.calculate_pcx_chord_emission(impact_factor, temp, w, mass, Lnu, v_outer, rmax=rmax, nr=500,
                                                         nlambda=2000, Lne=4.0, R_outer=r_outer)

    L = 150.0 / 0.004
    d = 0.88
    F = 20.7

    nx = 6000
    ny = 4000

    data = np.zeros((ny, nx))
    x = np.arange(1, nx+1, 1)
    y = np.arange(1, ny+1, 1)

    x0 = 3000
    y0 = 2000

    XX, YY = np.meshgrid(x, y)
    R = np.sqrt((XX-x0)**2 + (YY-y0)**2)

    nprocs = 4
    sub_partitions = 5000
    max_partitions = nprocs * sub_partitions

    split_r = np.array_split(R.flatten(), max_partitions)
    results = [0.0 for _ in split_r]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for subpartition in range(sub_partitions):
            print(subpartition)
            fut_to_spl = dict()
            for i in range(nprocs):
                idx = nprocs*subpartition + i
                # print(idx, subpartition, i)
                split = split_r[idx]
                future = executor.submit(models.general_model, split, L, d, F, warr, spectrum)
                fut_to_spl[future] = idx

            for future in concurrent.futures.as_completed(fut_to_spl, timeout=30):
                current_idx = fut_to_spl[future]
                output = future.result()
                # print(output.shape)
                results[current_idx] = output

    results = np.concatenate(results)
    results.shape = (ny, nx)

    fig, ax = plt.subplots()
    ax.imshow(results)
    plt.show()

if __name__ == "__main__":
    w0 = 487.98634
    mu = 39.948
    Ti = 0.5
    Vouter = 4000.0
    Router = 35.0
    Rmax = 42.0
    impact_fac = 33.0
    mom_dif_length = 20.0

    main(w0, mu, Ti, Vouter, Router, impact_fac, mom_dif_length, rmax=Rmax)
    
    