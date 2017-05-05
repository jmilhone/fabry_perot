from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import forward_model as fm
import time
import plottingtools.core as ptools
import cPickle as pickle
import pymultinest


rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

def eval_airy(wavelength, cos_theta, d, F):
    Q = 2. * F**2 / np.pi
    airy = 1.0 / (1.0 + Q * np.sin(np.pi * 2.e6 * d * cos_theta / wavelength)**2)
    return airy

def eval_spec(wavelength, amp, w0, sigma):
    norm = 1.0 / sigma / np.sqrt(2.0 * np.pi)
    exp = np.exp(-0.5*(wavelength-w0)**2 / sigma**2)
    return amp * norm * exp

def forward(r, L, d, F, Ti, mu, w0, nlambda=512):
    sigma = 3.276569e-5 * np.sqrt(Ti/mu) * w0
    w_arr = np.linspace(w0 - 5*sigma, w0 + 5*sigma, nlambda)
    nr = len(r)

    ll = np.tile(w_arr, (nr, 1)).T

    cos_th = L / np.sqrt(L**2 + r**2)
    cth = np.tile(cos_th, (len(w_arr), 1))


    spec = eval_spec(w_arr, 1.0, w0, sigma) 
    spec = np.tile(spec, (nr, 1)).T

    linear_out = np.trapz(spec*eval_airy(ll, cth, d, F), w_arr, axis=0)
    linear_out /= linear_out.max()
    return linear_out

def forward2(r, L, d, F, Ti, mu, amp, w0, nlambda=512):
    sigma = [3.276569e-5 * np.sqrt(Ti[idx]/mu[idx]) * w for idx, w in enumerate(w0)]
    minW = min(w0)
    maxW = max(w0)
    max_sigma = max(sigma)

    w_arr = np.linspace( minW- 5.0*max_sigma, maxW+5.0*max_sigma, nlambda)
    
    spec = 0.0
    for idx, w in enumerate(w0):
        spec += eval_spec(w_arr, amp[idx], w, sigma[idx])
    nr = len(r)

    ll = np.tile(w_arr, (nr, 1)).T

    cos_th = L / np.sqrt(L**2 + r**2)
    cth = np.tile(cos_th, (len(w_arr), 1))

    spec = np.tile(spec, (nr, 1)).T
    linear_out = np.trapz(spec*eval_airy(ll, cth, d, F), w_arr, axis=0)

    #plt.plot(r**2, linear_out)
    #plt.show()

    return linear_out

if __name__ == "__main__":
    #d =0.884157316074 
    #L = 150.535647184 / 0.004
    #Q = 20.0

    analyzer = pymultinest.Analyzer(n_params=8, outputfiles_basename="saves/full_solver_run4/fp_full_")
    stats = analyzer.get_mode_stats()
    mode = stats['modes'][0]
    
    mode_vals = mode['maximum a posterior']
    L = mode_vals[0]
    d = mode_vals[1]
    F = mode_vals[2]
    Ti0 = mode_vals[3]
    Ti1 = mode_vals[4]
    Amp_Th = mode_vals[5]
    Amp_Ar = mode_vals[6]
    rscale = mode_vals[7]
    #d = 8.82693201e-01
    #L = 1.50481176e+02/ .004
    #F = 1.92960567e+01 

    lambda_0 = 487.873302
    lambda_1 = 487.98634

    #Ti0 = 1000.0 * .025 / 300.0   #1000 K in eV
    #Ti0 = 1.10068828e+03 * .025 / 300.0
    mu0 = 232.0
    #Ti1 = 1000.0 * .025 / 300.0   #1000 K in eV
    #Ti1 = 5.50047357e-01
    mu1 = 40.0 
    
    
    #sigma0 = 3.276569e-5 * np.sqrt(Ti0/mu0) * lambda_0
    #sigma1 = 3.276569e-5 * np.sqrt(Ti1/mu1) * lambda_1


    #lambda0_arr = np.linspace(lambda_0 - 5*sigma0, lambda_0 + 5*sigma0, 1000)
    #lambda1_arr = np.linspace(lambda_1 - 5*sigma1, lambda_1 + 5*sigma1, 1000)

    r_arr = np.load("rbins2.npy")
    print r_arr
    print len(r_arr)
    print Amp_Ar / Amp_Th
    linear_out = forward2(r_arr, L, d, F, [Ti0, Ti1], [mu0, mu1], [Amp_Th , Amp_Ar], [lambda_0, lambda_1], nlambda=512)

    #up = [1012, 1287, 2981, 3260, 4957, 5238, 6935, 7218, 8921, 9202]
    #down = [1096, 1451, 3059, 3417, 5030, 5386, 7001, 7357, 8950, 9310]


    #up = [629.623347725, 714.720169017, 1085.42826571, 1137.60395569, 1400.52574414, 1442.16639817, 1656.98057925, 1693.0267334 ]
    #down = [665.805917667,759.449432155, 1106.4501977, 1165.20296944, 1417.16065427, 1462.94954117, 1671.30219889, 1709.60368507]

    #inds = []
    #for i in range(8):
    #    j = np.argmin(np.abs(r_arr-up[i]))
    #    k = np.argmin(np.abs(r_arr-down[i]))
    #    inds += range(j,k+1)
    ##r_arr = r_arr[inds]
    #nr = len(r_arr)

    #t0 = time.time()
    #linear_out0 = forward(r_arr, L, d, F, Ti0, mu0, lambda_0)
    #linear_out1 = forward(r_arr, L, d, F, Ti1, mu1, lambda_1)
    ##ll0 = np.tile(lambda0_arr, (nr, 1)).T
    ##ll1 = np.tile(lambda1_arr, (nr, 1)).T
    ##
    ##cos_th = L / np.sqrt(L**2 + r_arr**2)
    ##cth = np.tile(cos_th, (len(lambda0_arr), 1))
    ##
    ##
    ##spec0 = eval_spec(lambda0_arr, 1.0, lambda_0, sigma0) 
    ##spec0 = np.tile(spec0, (nr, 1)).T
    ##spec1 = eval_spec(lambda1_arr, 1.1, lambda_1, sigma1) 
    ##spec1 = np.tile(spec1, (nr, 1)).T
    ##
    ##
    ##linear_out0 = np.trapz(spec0*eval_airy(ll0, cth, d, Q), lambda0_arr, axis=0)
    ##linear_out1 = np.trapz(spec1*eval_airy(ll1, cth, d, Q), lambda1_arr, axis=0)

    ## mess with amplitudes here
    ##linear_out0 /= linear_out0.max()
    ##linear_out0 *= 0.5
    ##linear_out1 /= linear_out1.max()

    ##linear_out0 *= 0.5
    ##rscale = 3600.0
    ##rscale =4.43328262e+03
    ##linear_out0 *= 7.77008369e+06
    ##linear_out1 *= 1.26358903e+07
    #linear_out0 *= Amp_Th# * 1.27 / 1.48
    #linear_out1 *= Amp_Ar# * 1.27 / 1.48
    ##rscale = 4050.0
    fall_off = np.exp(-(r_arr / rscale)**2)
    ##plt.plot(r_arr, linear_out0*fall_off, 'r', label="{}".format(lambda_0))
    ##plt.plot(r_arr, linear_out1*fall_off, 'b', label="{}".format(lambda_1))
    #vals = fall_off*(linear_out0 + linear_out1)
    vals = fall_off*linear_out
    #t1 = time.time()
    #print t1-t0, "seconds"
    with open("saves/full_solver_run4/fp_ringsum_params.p") as infile:
        data = pickle.load(infile)
    #print data.keys()
    fig, ax = plt.subplots()
    ax.plot(data['binarr']**2, data['ringsum'], 'r')
    ax.plot(r_arr**2, vals)
    #ptools.add_thick_box(ax)
    ##plt.plot(r_arr, linear_out0, 'r')
    ##plt.plot(r_arr, linear_out1, 'b')
    ##lg = plt.legend()
    ##ax.set_xlim(0, 2000)
    ##ax.set_ylim(0, 1.0)
    #ptools.add_labels(ax, r"R${}^2$ (px${}^2$)", "Counts (A.U.)")
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    #plt.tight_layout()
    ##for i in range(8):
    ##    plt.axvline(x=up[i], color='r')
    ##    plt.axvline(x=down[i], color='g')
    #fig.savefig("forward_model_example.pdf")
    #plt.show()

    #darr = np.linspace(.88, .89, 100000)
    #m0 = 2.e6 * darr / lambda_0
    #r0 = np.sqrt( m0**2 / np.floor(m0)**2 -1.0)
    #
    #m1 = 2.e6 * darr / lambda_1
    #r1 = np.sqrt( m1**2 / np.floor(m1)**2 -1.0)
    #
    #plt.plot(darr, r0, label="{}".format(lambda_0))
    #plt.plot(darr, r1, label="{}".format(lambda_1))
    #plt.legend()
    #plt.show()
    #ll0 = np.tile(lambda_arr0, (len(r_arr), 1)).T
    #t0 = time.time()
    #cos_th = L / np.sqrt(L**2 + r_arr**2)
    #
    ##ll, cth = np.meshgrid(lambda_arr, cos_th, indexing='ij')
    #
    #cth = np.tile(cos_th, (len(ll), 1))
    #
    #
    #print cth.shape
    #print ll.shape
    #Ti = 1000.0 * .025 / 300.0   #1000 K in eV
    #mu = 232.0
    #sigma = 3.276569e-5 * np.sqrt(Ti/mu) * lambda_0
    #
    #spec = eval_spec(lambda_arr, 1.0, lambda_0, sigma) 
    #
    #spec = np.tile(spec, (len(cos_th), 1)).T
    #
    #linear_out = np.trapz(spec*eval_airy(ll, cth, d, 20.0), lambda_arr, axis=0)
    #
    #
    ##half = np.isclose(linear_out, 0.5, atol=0.005)
    ##print half
    ##inds = [i for i, x in enumerate(half) if x]
    #print len(inds), len(linear_out)
    #fall_off = np.exp(-r_arr / (5*r_arr.max()))
    #print time.time()-t0, "seconds"
    #
    #ncth = cos_th[inds[0]:inds[1]]
    #print inds[0], inds[1]
    #npts = len(ncth)
    #nll = np.tile(lambda_arr, (npts, 1)).T
    #ncth = np.tile(ncth, (len(lambda_arr), 1))
    #
    #neval = eval_spec(lambda_arr, 1.0, lambda_0, sigma) 
    #nspec = np.tile(neval, (npts,1)).T
    #print nspec.shape
    #print nll.shape, ncth.shape
    #print eval_airy(nll, ncth, d, 20.0).shape
    #nlinear_out = np.trapz(nspec*eval_airy(nll, ncth, d, 20.0), lambda_arr, axis=0)
    ##plt.plot(r_arr, linear_out*fall_off)
    #plt.plot(r_arr, linear_out)
    #plt.plot(r_arr[inds[0]:inds[1]], nlinear_out, 'g')
    plt.show()





























