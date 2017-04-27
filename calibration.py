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

np.seterr(invalid='ignore')  # I got sick of invalid values that happening during minimization
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
global_L_limits = (37000, 38000)
global_d_limits = (.87, .89)

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
def ringsum2(R, weights, m0, L, d, peaks, dlambda_arr, out=None, index=None):

    #ringsums = np.zeros( (len(peaks), len(dlambda_arr)-1) )
    ringsums = np.zeros( (1, len(dlambda_arr)-1) )
    for j, peak in enumerate(peaks[0:1]):
        sq = np.sqrt(1.0 + (peak/L)**2)
        num = 2.e6 * d * sq
        denom = 2.e6 * d + dlambda_arr * (m0 - j) * sq
        temp = num / denom
        rarr = L * np.sqrt(temp**2 - 1.0)

        # rarr is the right direction, but the histogram requires the bins increase
        rmin = np.nanmin(rarr)
        rmax = np.nanmax(rarr)
        #t0 = time.time()
        inds = np.where(np.logical_and(R >= rmin, R< rmax))
        #print time.time()-t0
        #sig, _ = np.histogram(R[inds], bins=rarr[::-1], weights=weights[inds])
        sig, _ = np.histogram(R[inds], bins=rarr, weights=weights[inds])
        sig = sig[::-1]  # correct the order again
        ringsums[j, :] = sig

    if out:
        print "exiting..."
        out.put((index, ringsums))
    else:
        return ringsums



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
    data = im.get_image_data(f, bgname=f_bg, color='b')
    # im.quick_plot(data)
    # plt.show()

    times += [time.time()]
    print "Image done reading, {0} seconds".format(times[-1] - times[-2])
    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)
    times += [time.time()]
    print "Center found, {0} seconds".format(times[-1] - times[-2])

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    new_binarr = np.concatenate(([0.0], binarr))
    n = len(new_binarr)

    rarr = binarr



    if L is None or d is None:
        peaks = fitting.peak_and_fit2(rarr, sigarr, thres=0.3, plotit=False, smooth_points=10)
        th_peaks = peaks[0::2]
        ar_peaks = peaks[1::2]

        print "Ar Peak locations: ", ar_peaks
        print "Th Peak locations: ", th_peaks

        analysis = mns.solve_L_d(save_dir, th_peaks, ar_peaks, L_lim=(149.0/.004, 151.0/.004) )
        L, d, Lsd, dsd = mns.L_d_results(analysis, th_peaks, ar_peaks, plotit=False)

        # Temporary checking V on Ar plasma
        #newf = "Images/0015676_000.nef"
        #newfbg = "Images/0015676_001.nef"
        #newf = "Images/0015406_000.nef"
        #newfbg = "Images/0015406_001.nef"

        #newf = "Images/0010105_000.nef"
        #newfbg = "Images/0010105_001.nef"
        #newdata = im.get_image_data(newf, newfbg, color='b') 
        #x1, y1 = rs.locate_center(newdata, x0, y0, binsize=binsize, plotit=False)
        #nbinarr, nsigarr = rs.quick_ringsum(newdata, x1, y1, binsize=binsize, quadrants=False)

        post = analysis.get_equal_weighted_posterior()
        Lpost = post[:, 0]
        dpost = post[:, 1]
        mpost = 2.e6 * dpost / c_lambda

        #arpeaks = fitting.peak_and_fit2(nbinarr, nsigarr, thres=.3, plotit=True, smooth_points=10)

        #pk0 = ar_peaks[0]
        dlambda_arr = np.arange(-n_dlambda/2, n_dlambda/2 + 1, 1) * delta_lambda
        #for j, pk0 in enumerate(ar_peaks):
        #    r = np.zeros_like(lam_arr)
        #    rstd = np.zeros_like(lam_arr)
        #    for idx, lam in enumerate(lam_arr):
        #        temp = 2.e6 * dpost * np.sqrt(Lpost**2 + pk0**2) / (2.e6*dpost*Lpost + lam * (np.floor(mpost) - j) * np.sqrt(Lpost**2 + pk0**2))
        #        rtemp = L * np.sqrt((temp)**2 -1.0)

        #        r[idx] = np.mean(rtemp)
        #        rstd[idx] = np.std(rtemp)
        #    plt.errorbar(lam_arr, r, yerr=rstd)
        #    #plt.plot(lam_arr, rstd / r * 100.0)
        #    plt.title("{0:d}".format(j))
        #    plt.show()
        #for j, pk0 in enumerate(ar_peaks):
        #    vpost = 2.998e8*(1.0 - mpost * Lpost / (np.floor(mpost)-1.0*j) / np.sqrt(Lpost**2 + pk0**2))
        #    rj = Lpost * np.sqrt( mpost**2 / (np.floor(mpost) - j)**2 -1 )
        #    plt.hist(rj, bins=75, normed=True)
        #    #plt.title("Peak {0:d}, V = {1:3.1f} +/- {2:3.2f} m/s".format(j,np.mean(vpost), np.std(vpost)))
        #    plt.title("Peak {0:d}, rj = {1:f} +/- {2:f} px".format(j, np.mean(rj), np.std(rj)))
        #    plt.show()

        ## Now run with finer resolution
        #save_dir = save_dir + "_refined"
        #analysis = mns.solve_L_d(save_dir, th_peaks, ar_peaks, L_lim=(L-5.0*Lsd, L+5.0*Lsd), d_lim=(d-5.0*dsd, d+5.0*dsd))
        #L, d, Lsd, dsd = mns.L_d_results(analysis, th_peaks, ar_peaks)
        #save_dir = save_dir + "_1"
        #analysis = mns.solve_L_d(save_dir, th_peaks, ar_peaks, L_lim=(L-1.0*Lsd, L+1.0*Lsd), d_lim=(d-0.0001, d+.0001))
        #L, d, Lsd, dsd = mns.L_d_results(analysis, th_peaks, ar_peaks)

    if gas == 'He':
        peaks = fitting.peak_and_fit2(rarr, sigarr, thres=.65, plotit=True, smooth_points=2)
    elif gas == 'Ar':
        peaks = fitting.peak_and_fit2(rarr, sigarr, thres=0.3, plotit=False, smooth_points=10)
        print peaks
        peaks = peaks[0::2]
        print peaks
    else:
        print "Can't analyze this gas"
        return None

    m0 = 2.0e6 * d / c_lambda
    ny, nx = data.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)
    xx, yy = np.meshgrid(1. * x - x0, 1. * y - y0)
    R = np.sqrt(xx ** 2 + yy ** 2)

    # Flatten data for simplicity
    R = R.flatten()
    data = data.flatten()
    isort = np.argsort(R)
    R = R[isort]
    data = data[isort]

    all_ringsums = []
    npts = len(Lpost)
    rv_dlambda_arr = dlambda_arr[::-1]

    nproc = 16 
    indexes = range(0, npts, 10)
    allringsums = np.zeros((len(th_peaks[0:4]), len(dlambda_arr)-1, len(indexes)))
    n, m = divmod(len(indexes), nproc)
    print n, m, len(indexes) * 1.0 / n
    t00 = time.time()
    for i in range(n):
        t0 = time.time()
        procs = []
        out = mp.Queue()
        if i != n-1:
            #print "in the initial if"
            for index in range(i*nproc, (i+1)*nproc):
                post_index = indexes[index]
                ring_index = index# + i*nproc
                #print ring_index, post_index, mpost[post_index], Lpost[post_index], dpost[post_index]
                p = mp.Process(target=rs.ringsum2, 
                        args=(R, data, np.floor(mpost[post_index]), Lpost[post_index], dpost[post_index], th_peaks[0:4], rv_dlambda_arr),
                        kwargs={'out': out, 'index': ring_index})
                procs.append(p)
                p.start()


            for j in range(nproc):
                #print "I'm here"
                tup = out.get()
                #print tup[0], tup[1].shape
                #print allringsums.shape, tup[0], tup[1][0,:].shape
                allringsums[:, :, tup[0]] += tup[1]
                #print "done adding ring sum"

            for p in procs:
                #print "Now im in the join loop"
                p.join()
                #print "done joining"
        else:
            for index in range(i*nproc, len(indexes)):
                #print index, len(indexes), i*nproc
                #print indexes[index]
                post_index = indexes[index]
                ring_index = index# + i*nproc
                #print ring_index, post_index, mpost[post_index], Lpost[post_index], dpost[post_index]
                p = mp.Process(target=rs.ringsum2, 
                        args=(R, data, np.floor(mpost[post_index]), Lpost[post_index], dpost[post_index], th_peaks[0:4], rv_dlambda_arr),
                        kwargs={'out': out, 'index': ring_index})
                procs.append(p)
                p.start()

            for j in range(len(procs)):
                #print "I'm here"
                tup = out.get()
                #print j, tup[0], tup[1].shape
                allringsums[:, :, tup[0]] += tup[1]
                #print "we got the data"

            for p in procs:
                #print "Now im in the join loop"
                p.join(5.0)
                if p.is_alive():
                    p.terminate()
                #print "after join statement" 

        print i, time.time() - t0
    print time.time() - t00, "seconds"

    for j in range(4):
        fig, ax = plt.subplots()
        for i in range(len(indexes)):
            ax.plot(dlambda_arr[0:-1],allringsums[j, :, i])
        ax.set_title("Order {}".format(j))
    plt.show()

    nrow, ncol, npost = allringsums.shape
    hwhm = np.zeros((nrow, npost))
    offset = np.zeros((nrow, npost)) 

    iL = np.abs((dlambda_arr + .006)).argmin()
    iR = np.abs((dlambda_arr - .006)).argmin()
    ind = range(iL, iR)
    i = np.abs(dlambda_arr).argmin()

    Ti = 1000.0 * .025 / 300.0
    sigma = 3.265e-5 * c_lambda * np.sqrt(Ti / lampmu)

    for order in range(nrow):
        for idx in range(npost):
            rsum = allringsums[order, :, idx]
            thres = .5 * rsum[i]
            try:
                iL = np.max(np.where(rsum[0:i] < thres))
                iR = np.min(np.where(rsum[i:] < thres)) + i
                inds = range(iL, iR+1)

                pk = np.trapz(dlambda_arr[inds] * rsum[inds], x=dlambda_arr[inds]) / np.trapz(rsum[inds], x=dlambda_arr[inds])

                maxval = np.max(rsum[ind])
                xdata = dlambda_arr[ind]
                ydata = rsum[ind] / maxval
                offset_guess = np.mean(rsum[-10:]) / maxval
                print offset_guess
                #opt, cov = curve_fit(lambda x, a, b, c: fitting.voigt_profile(x, a, b, sigma, pk)+c,
                #                     xdata, ydata, p0=[0.01, 0.003, offset_guess])
                opt, cov = curve_fit(lambda x, a, b: fitting.voigt_profile(x, a, b, sigma, pk),
                                     xdata, ydata, p0=[0.01, 0.003])
                #print opt[1]
                hwhm[order, idx] = opt[1]
                #offset[order, idx] = opt[2] * maxval
            except ValueError, e:
                print "There is something wrong with order {}, post {}".format(order, idx)
                print e
    inds = np.where(hwhm[0, :] > 0.0) 

    plt.figure()
    plt.hist(np.squeeze(hwhm[0, inds]), bins=20)
    plt.show()
    #plt.figure()
    #plt.hist(np.squeeze(offset[0,inds]), bins=20)
    #plt.show()
    #for i in range(0, npts, 5):
    #    #t0 = time.time()
    #    all_ringsums.append(rs.ringsum2(R, data, np.floor(mpost[i]), Lpost[i], dpost[i], th_peaks, rv_dlambda_arr))
    #    #print all_ringsums[-1].shape, dlambda_arr.shape, rv_dlambda_arr.shape
    #    #print i, time.time() - t0 
    ###Find fitting region:
    #iL = np.abs((dlambda_arr + .006)).argmin()
    #iR = np.abs((dlambda_arr - .006)).argmin()
    #ind = range(iL, iR)

    #Ti = 1000.0 * .025 / 300.0
    #sigma = 3.265e-5 * c_lambda * np.sqrt(Ti / lampmu)

    #hwhm = []

    #for rsum_ in all_ringsums:
    #    rsum = rsum_[0, :]
    #    #rsum = all_ringsums[i][0]
    #    i = np.abs(dlambda_arr).argmin()
    #    thres = .5 * rsum[i]
    #    iL = np.max(np.where(rsum[0:i] < thres))
    #    iR = np.min(np.where(rsum[i:] < thres)) + i
    #    inds = range(iL, iR+1)
    #    pk = np.trapz(dlambda_arr[inds] * rsum[inds], x=dlambda_arr[inds]) / np.trapz(rsum[inds], x=dlambda_arr[inds])
    #    #print pk

    #    maxval = np.max(rsum[ind])
    #    xdata = dlambda_arr[ind]
    #    ydata = rsum[ind] / maxval

    #    opt, cov = curve_fit(lambda x, a, b: fitting.voigt_profile(x, a, b, sigma, pk),
    #                         xdata, ydata, p0=[0.01, 0.003])
    #    print "HWHM gamma= {0} nm".format(opt[1])
    #    hwhm.append(opt[1]) 

    #plt.hist(hwhm, normed=True, bins=20)
    #plt.show()


        #plt.plot(dlambda_arr[0:-1], rsum[0,:])
    #plt.show()
    #ringsums, lambda_arr = rs.ringsum(R, data, m0, L, d, peaks, delta_lambda, ndl=n_dlambda)

    #w_peaks = []
    #for rsum in ringsums:
    #    i = np.abs(lambda_arr).argmin() # close enough to the peak
    #    thres = .5 * rsum[i]
    #    iL = np.max(np.where(rsum[0:i] < thres))
    #    iR = np.min(np.where(rsum[i:] < thres)) + i
    #    inds = range(iL, iR+1)
    #    pk = np.trapz(lambda_arr[inds] * rsum[inds], x=lambda_arr[inds]) / np.trapz(rsum[inds], x=lambda_arr[inds])
    #    w_peaks += [pk]


    ###Find fitting region:
    #iL = np.abs((lambda_arr + .006)).argmin()
    #iR = np.abs((lambda_arr - .006)).argmin()
    #ind = range(iL, iR)

    #Ti = 1000.0 * .025 / 300.0
    #sigma = 3.265e-5 * c_lambda * np.sqrt(Ti / lampmu)
    #
    #hwhm = []
    #for i, pk in enumerate(w_peaks):
    #    maxval = np.max(ringsums[i][ind])
    #    xdata = lambda_arr[ind]
    #    ydata = ringsums[i][ind] / maxval

    #    opt, cov = curve_fit(lambda x, a, b: fitting.voigt_profile(x, a, b, sigma, pk),
    #                         xdata, ydata, p0=[0.01, 0.0036])
    #    print "Peak {0}: HWHM gamma= {1} nm".format(i, opt[1])
    #    hwhm.append(opt[1]) 

    #    voigt = fitting.voigt_profile(lambda_arr, opt[0], opt[1], sigma, pk)

    #    fig, ax = plt.subplots()
    #    ax.plot(lambda_arr, maxval*voigt, 'r')
    #    ax.plot(lambda_arr, ringsums[i], 'b')
    #    plt.show()

    #return L, d, hwhm, x0, y0

if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    save_dir = "June_09_2016_Calibration_test1"

    # Old calibration used to flow
    #fname = join(folder, "Th_lamp_488_calib_5m.nef")
    #bg_fname = join(folder, "Th_lamp_488_calib_5m_background.nef")
    #save_dir = "old_flow_shots"

    center_guess = (3068.56, 2033.17)
    #L, d, hwhm, x0, y0 = run_calibration(fname, bg_fname, center_guess, save_dir, gas='Ar')
    run_calibration(fname, bg_fname, center_guess, save_dir, gas='Ar')

    #fname = join(folder, "thorium_5_min_he_3.nef")
    #bg_fname = join(folder, "thorium_5_min_he_bg.nef")
    # center_guess = (3068.39, 2031.85)
    # He guess
    #center_guess = (3040.05627213, 2024.06787634)
   # L, d, hwhm, x0, y0 = run_calibration(fname, bg_fname, center_guess, save_dir, gas='Ar')

