from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import minimize_scalar, fmin, curve_fit, differential_evolution
import time
import cPickle as pickle
from analysis.datahelpers import smooth
import multinest_solver as mns
import multiprocessing as mp
import model
from scipy.stats import norm

#np.seterr(invalid='ignore')  # I got sick of invalid values that happening during minimization
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

Ar_params = {
    'binsize': 0.1,
    'c_lambda': 487.873302,
    'delta_lambda': 0.0001*4,
    'duc': 0.88,
    'n_dlambda': 512 /4 ,
    'Tlamp': 1000.0 * .025 / 300.0,  # .025 eV per 300 K
    'lampmu': 232.0,
}
He_params = {
    'binsize': 0.1*4.0,
    'c_lambda': 468.6195,
    'delta_lambda': 0.0001*4.0,
    'duc': 0.88,
    'n_dlambda': 512/4,
    'Tlamp': 1000.0 * .025 / 300.0,  # .025 eV per 300 K
    'lampmu': 232.0,
}

def run_calibration(f, f_bg, center_guess, save_dir, gas='Ar', L=None, d=None):
    times = [time.time()]
    x0, y0 = center_guess

    if gas == 'Ar':
        binsize = Ar_params['binsize']
        c_lambda = Ar_params['c_lambda']
        delta_lambda = Ar_params['delta_lambda']
        duc = Ar_params['duc']
        n_dlambda = Ar_params['n_dlambda']
        Tlamp = Ar_params['Tlamp']
        lampmu = Ar_params['lampmu']
    elif gas == 'He':
        binsize = He_params['binsize']
        c_lambda = He_params['c_lambda']
        delta_lambda = He_params['delta_lambda']
        duc = He_params['duc']
        n_dlambda = He_params['n_dlambda']
        Tlamp = He_params['Tlamp']
        lampmu = He_params['lampmu']
    else:
        print "Did not recognize gas.  Exiting..."
        return None
    # f_bg=None
    print f[-3:].lower()
    if f[-3:].lower() == "nef":
        #print "changed it to green dumbass"
        #data = im.get_image_data(f, bgname=f_bg, color='g')
        data = im.get_image_data(f, bgname=f_bg, color='b')
    elif f[-3:].lower() == "npy":
        data = np.load(f)
    else:
        print "need a valid file"
        return
    im.quick_plot(data)
    plt.show()

    times += [time.time()]
    print f
    print "Image done reading, {0} seconds".format(times[-1] - times[-2])
    x0, y0 = rs.locate_center(data, x0, y0, binsize=0.1, plotit=False)
    fig, ax = plt.subplots(2)
    ax[0].plot(data[int(np.ceil(y0)), :])
    ax[1].plot(data[int(np.floor(y0)), :])
    plt.show()
    times += [time.time()]
    print "Center found, {0} seconds".format(times[-1] - times[-2])
    
    binarr, UL, UR, BL, BR = rs.quick_ringsum(data, x0, y0, binsize=1.0)


    r = np.concatenate(([0.0], binarr))
    rr = 0.5 * (r[0:-1] + r[1:])
    rdiff = np.diff(r)
    fig, ax = plt.subplots()
    plt.plot(binarr**2, UL, label='UL')
    plt.plot(binarr**2, UR, label='UR')
    plt.plot(binarr**2, BL, label='BL')
    plt.plot(binarr**2, BR, label='BR')
    #plt.plot(binarr**2, UL/UL.max(), label='UL')
    #plt.plot(binarr**2, UR/UR.max(), label='UR')
    #plt.plot(binarr**2, BL/BL.max(), label='BL')
    #plt.plot(binarr**2, BR/BR.max(), label='BR')
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax.legend(loc='upper right')
    plt.show()
    

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=1.0, quadrants=False)
    bg_Th = np.abs(norm(scale=0.2).rvs(data.shape))
    _, bg_sig_Th = rs.quick_ringsum(bg_Th, x0, y0, binsize=1.0, quadrants=False)

    sigarr -= bg_sig_Th

    _, err_sqd = rs.quick_ringsum(100.0**2 * np.ones_like(data), x0, y0, binsize=1.0, quadrants=False)
    err = np.sqrt(err_sqd)
    fig, ax = plt.subplots()
    ax.errorbar(rr, sigarr, xerr=rdiff/2.0, yerr=err)
    rrr = np.zeros_like(rr)
    for idx, rval in enumerate(rr):
        rrr[idx] = norm(scale=rdiff[idx]/6.0, loc=rval).rvs(1)[0]
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax.plot(rrr, sigarr,'g')
    plt.show()

    binarr = np.concatenate(([0.0], binarr))
    binarr = 0.5*(binarr[0:-1] + binarr[1:])
    #ny, nx = data.shape
    #x = np.arange(0, nx, 1)
    #y = np.arange(0, ny, 1)
    #xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    #RR = np.sqrt(xx**2 + yy**2)

    #rr = np.linspace(0, 2000.0, 2000)
    #tot_areas = np.pi * rr**2
    #areas = np.diff(tot_areas)

    #ss, _ = np.histogram(RR, bins=rr, weights=data)
    #print rr.shape
    #print ss.shape
    #rr = rr[1:]
    #ss = ss / areas
    #fig, ax = plt.subplots()
    #ax.plot(binarr, sigarr, 'g')
    #ax1 = ax.twinx()
    #ax1.plot(rr, ss, 'b')
    #plt.show()

    #print "Subtracting min value from ring sum"
    #sigarr -= sigarr.min()

    print len(binarr)
    #np.save("rbins2.npy",binarr)
    #np.save("old_thorium_data.npy", data)
    new_binarr = np.concatenate(([0.0], binarr))
    LL, RR = fitting.get_peaks(binarr, sigarr, thres1=0.2,thres2=0.3, npks=8)
    z = fitting.peak_and_fit2(binarr, sigarr, thres=0.2, plotit=False)
    ampsAr = []
    ampsTh = []
    locAr = []
    locTh = []
    for i in range(3):
        temp1 = np.max(sigarr[LL[2*i]:RR[2*i]])
        temp2 = np.max(sigarr[LL[2*i+1]:RR[2*i+1]])
        ampsAr.append(temp2)
        ampsTh.append(temp1)
        locAr.append(z[2*i+1])
        locTh.append(z[2*i])

    print "sig max",sigarr.max()
    ampsAr = np.array(ampsAr)
    ampsTh = np.array(ampsTh)
    locAr = np.array(locAr)
    locTh = np.array(locTh)

    print ampsAr
    print ampsTh
    relamps = [ampsAr[i]/ampsTh[i] for i in range(3)]
    relA = np.mean(relamps)
    print relamps, relA
    print z
   
    amps = np.hstack((ampsAr, relA*ampsTh))
    locs = np.hstack((locAr, locTh))
    pfit = np.polyfit(locs**2, np.log(amps),1)
    print "r0 fit: ",1.0 / np.sqrt(np.abs(pfit[0]))
    r0_fit =  1.0 / np.sqrt(np.abs(pfit[0]))
    print "sig max",sigarr.max()
    plt.plot(locTh**2, ampsTh, 'o')
    plt.plot(locTh**2, relA*ampsTh, 'o')
    plt.plot(locAr**2, ampsAr, 'o')
    plt.plot(binarr**2, sigarr)
    plt.plot(binarr**2, np.exp(np.polyval(pfit,binarr**2)))
    plt.show()


    #for i in range(8):
    #    print L[i], R[i]
    inds = [] 
    print "Changed peaks from 8 to 6"
    print "Using all indexes between argon and th peaks"
    for i in range(3):
        inds += range(LL[2*i], RR[2*i+1])
       
    amp0 = np.max(sigarr[LL[0]:RR[0]])
    amp1 = np.max(sigarr[LL[1]:RR[1]])

    i = range(LL[0], RR[0])
    pk0 = np.trapz(binarr[i] * sigarr[i], x=binarr[i]) / np.trapz(sigarr[i], x=binarr[i])

    i = range(LL[1],RR[1])
    pk1 = np.trapz(binarr[i] * sigarr[i], x=binarr[i]) / np.trapz(sigarr[i], x=binarr[i])
    print pk0, pk1
    r = binarr[inds]

    d =0.884157316074 
    L = 150.535647184 / 0.004
    Q = 20.0

    lambda_0 = 487.873302
    lambda_1 = 487.98634

    Ti0 = 1000.0 * .025 / 300.0   #1000 K in eV
    mu0 = 232.0
    Ti1 = 1000.0 * .025 / 300.0   #1000 K in eV
    mu1 = 40.0 

    sigma0 = 3.276569e-5 * np.sqrt(Ti0/mu0) * lambda_0
    sigma1 = 3.276569e-5 * np.sqrt(Ti1/mu1) * lambda_1
    linear_out = model.forward2(r, L, d, Q, [Ti0, Ti1], [mu0, mu1], [amp0, amp1], [lambda_0, lambda_1])
    #linear_out0 = model.forward(r, L, d, Q, Ti0, mu0, lambda_0)
    #linear_out1 = model.forward(r, L, d, Q, Ti1, mu1, lambda_1)
    #linear_out0 *= 0.5
    #linear_out0 *= 1.5e7
    #linear_out1 *= 1.5e7

    rscale = 4500.0
    fall_off = np.exp(-(r / rscale)**2)
    #linear_out0 *= amp0 / np.exp(-pk0/rscale)
    #linear_out1 *= amp1 / np.exp(-pk1/rscale)

    fig, ax = plt.subplots()
    plt.plot(binarr, sigarr, 'b')
    print len(binarr)
    #ax1 = ax.twinx()
    #ax1.plot(r, fall_off * (linear_out0 + linear_out1), 'r')
    #ax.plot(r, fall_off * (linear_out0 + linear_out1), 'r')
    ax.plot(r, fall_off * linear_out, 'r')
    for i in range(3):
        plt.axvline(x=binarr[LL[2*i]], color='g')
        plt.axvline(x=binarr[RR[2*i+1]], color='r')
    plt.show()


    finesse_range = [15, 25]
    #Ti_Ar_range = [300.0, 11000]
    Ti_Ar_range = [300.0, 11000/2.0]
    Ti_Th_range = [300.0, 3000]
    d_range = [.87, .89]
    L_range = [145.0, 155.0]
    a0 = amp0 / np.exp(-(pk0/rscale)**2)
    a1 = amp1 / np.exp(-(pk1/rscale)**2)
    #Th_amp_range = [.5*a0, 2.0*a0]
    #Ar_amp_range = [.5*a1, 2.0*a1]
    Th_amp_range = [.75*a0, 1.5*a0]
    Ar_amp_range = [.75*a1, 1.5*a1]
    rscale_range = [500.0, 3500.0]
    print a0, a1 

    #print len(binarr)
    #plt.plot(binarr, sigarr, 'r')

    #theta_min = 0.0
    #L = 150.508755008
    #rmax = 1800.0
    #cos_theta_min = L / np.sqrt(L**2 + rmax**2)
    #cos_theta = np.linspace(cos_theta_min,1.0, 1000.0)
    #newr = L * np.sqrt( 1.0 /cos_theta**2 -1.0)
    #newr = newr[::-1]
    #ny, nx = data.shape
    #x = np.arange(0, nx, 1)
    #y = np.arange(0, ny, 1)
    #print x0, y0
    #xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    #R = np.sqrt(xx**2 + yy**2)


    #newsig, _ = np.histogram(R.flatten(), bins=newr, weights=data.flatten())
    ##plt.plot(newr)
    #plt.plot(newr[1:], newsig,'b')
    #plt.show()
    #plt.plot(np.diff(newr),'.')
    #plt.show()
    #plt.plot(R.flatten(), data.flatten(), '.')
    #plt.show()
    #params = {"finesse_range": finesse_range,
    #          "Ti_Ar_range": Ti_Ar_range,
    #          "Ti_Th_range": Ti_Th_range,
    #          "d_range": d_range,
    #          "L_range": L_range,
    #          "Th_amp_range": Th_amp_range,
    #          "Ar_amp_range": Ar_amp_range,
    #          "rscale_range": rscale_range,
    #          "binarr": binarr,
    #          "ringsum": sigarr,
    #          "ringsum_sd": 100.0 + 0.01 * sigarr,
    #          "ind": inds}

    params = {"finesse_range": finesse_range,
              "Ti_Ar_range": Ti_Ar_range,
              "Ti_Th_range": Ti_Th_range,
              "d_range": d_range,
              "L_range": L_range,
              "Th_amp_range": Th_amp_range,
              "Ar_amp_range": Ar_amp_range,
              "rscale_range": rscale_range,
              "binarr": rr,
              "binsize": rdiff,
              "ringsum": sigarr,
              "ringsum_sd": 100.0 + 0.01 * sigarr,
              "ind": inds}

    #params = {"finesse_range": finesse_range,
    #          "Ti_Ar_range": Ti_Ar_range,
    #          "Ti_Th": 1000.0 * .025 / 300.0,
    #          "d_range": d_range,
    #          "L_range": L_range,
    #          "rel_amp": relA,
    #          "Ar_amp_range": Ar_amp_range,
    #          "rscale": r0_fit,
    #          "binarr": binarr,
    #          "ringsum": sigarr,
    #          "ringsum_sd": 100.0 + 0.03 * sigarr,
    #          "ind": inds}


    folder = mns.fix_save_directory(save_dir)
    with open(folder + "fp_ringsum_params.p", 'wb') as outfile:
        pickle.dump(params, outfile)


if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    #fname = join(folder, "thorium_ar_5_min_1.nef")
    #bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")

    #folder = "FM_saves"
    #fname = "model1.npy"
    #fname = "model1_0.02_noise.npy"
    #fname = join(folder,fname)
    #fname = "thorium_lamp_data.npy"
    #fname = "old_thorium_data.npy"
    #folder = "Images"
    #fname = join(folder, "Th_lamp_488_calib_5m.nef")
    #fname = join(folder, "0012443_000.nef")
    #bg_fname = join(folder, "Th_lamp_488_calib_5m_background.nef")
    #fname = join(folder, "thorium_ar_5_min_1.nef")
    #bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    #fname = join(folder, "DSC_0007.nef")
    bg_fname = None
    #bg_fname = join(folder, "0012443_001.nef")
    #fname = "thorium_lamp_data.npy"
    fname = join("synthetic_data/test1/","th_lamp_ar.npy") 
    #fname = join(folder, "DSC_0098.NEF")
    #fname = join(folder, "0025203_000.nef")
    #fname = "FM_150_88_195_4000_0.3eV.npy"
    #bg_fname = None
    #save_dir = "full_solver_run17"
    #save_dir = "new_full_solver_run0"

    #save_dir = "solver3_run16"
    #save_dir = "finesse_solver3"
    #center_guess = (3068.56, 2033.17)
    #center_guess = (3068.57005213, 2032.17646934)
    #center_guess = (3040.78197222, 2003.29777455) # Old thorium calibration data center
    center_guess = (3018.0, 2010.0)
    #center_guess = (3041., 2028.)
    #center_guess = (3068.56, 2033.17)
    #center_guess = (3040.78197222, 2003.29777455) # old Th calib center guess
    #center_guess = (3068.57005213, 2032.17646934)
    #center_guess = (3018.0, 2010.0)
    #center_guess = (3070.0, 2270.0)
    #center_guess = (3046, 2207)
    #center_guess = (3070.0, 2500.0)
    #folder ="ForwardModelData"
    #save_dir = "NewCalibration_Ar_run0"
    save_dir = "finesse_solver10"
    #fname = join(folder, "th_lamp_FM_14987_8824.npy")
    #bg_fname = None
    #save_dir = "FW_test_14987_8824_0"
    #center_guess = (3018, 2010)
    # Old calibration used to flow
    #fname = join(folder, "Th_lamp_488_calib_5m.nef")
    #bg_fname = join(folder, "Th_lamp_488_calib_5m_background.nef")
    #save_dir = "old_flow_shots"

    #L, d, hwhm, x0, y0 = run_calibration(fname, bg_fname, center_guess, save_dir, gas='Ar')
    run_calibration(fname, bg_fname, center_guess, save_dir, gas='Ar')

    #fname = join(folder, "thorium_5_min_he_3.nef")
    #bg_fname = join(folder, "thorium_5_min_he_bg.nef")
    # center_guess = (3068.39, 2031.85)
    # He guess
    #center_guess = (3040.05627213, 2024.06787634)
   # L, d, hwhm, x0, y0 = run_calibration(fname, bg_fname, center_guess, save_dir, gas='Ar')

