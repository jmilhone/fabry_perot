import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import fmin
import sys
sys.path.append('../')
from core.synthetic_data import full_pattern, forward_model
from core.ringsum import smAng_ringsum, locate_center
from tools.file_io import h5_2_dict
from tools.plotting import center_plot, my_hist, my_hist2d

def make_synthetic_data():
    full_pattern(204.38721, 3.875373, 4., 487.98634, 39.948, 0.15, 0.0, nprocs=30,plotit=True,saveit='./Data/synthetic/error_test_01.h5')

def ringsum(data, x0, y0, binsize=1.0):
    '''
    performs a ringsum given a center and data with 
    bins centered on pixel

    Args:
        data (np.ndarray): image data
        x0 (float): x center
        y0 (float): y center
        binsize (float, default=1.0): binsize
    Returns:
        binarr (np.ndarray): array of bins
        sig (np.ndarray): ringsum values
    '''
    binarr, sig = smAng_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    binarr = np.concatenate(([0.0], binarr))
    binarr = 0.5 * (binarr[0:-1] + binarr[1:])
    return binarr, sig

def perf_center():
    data = h5_2_dict('./Data/synthetic/error_test_01.h5')
    A = data['2Ddata']
    xg = []
    yg = []
    x0 = []
    y0 = []
    X,Y = np.meshgrid(np.arange(A.shape[1]),np.arange(A.shape[0]))
    ### 4200 takes 15.3 hours ###
    for i in range(4200):
        print '{0} of 4200'.format(i+1)
        xxg = np.random.rand() * 10. + 3008.-5.
        yyg = np.random.rand() * 10. + 2008.-5.
        A = np.copy(data['2Ddata'])
        A += np.abs(np.random.normal(scale=0.05, size=A.shape))
        xo = np.random.rand() * 2000. + 2008.
        yo = np.random.rand() * 2000. + 1008.
        sig = np.random.rand() * (2000. - 500.) + 500.
        amp = np.random.rand() * (3. - 1) + 1.
        off = amp * np.exp(-0.5*(X-xo)**2/sig**2)*np.exp(-0.5*(Y-yo)**2/sig**2)
        A += A*off
        fade = np.exp(-0.5*(X-3007.5)**2/1500**2)*np.exp(-0.5*(Y-2007.5)**2/1500**2)
        A += A*fade
        x, y = locate_center(A, xguess=xxg , yguess=yyg)
        x0.append(x)
        y0.append(y)
        xg.append(xxg)
        yg.append(yyg)
        del A

    x0 = np.array(x0)
    y0 = np.array(y0)
    xg = np.array(xg)
    yg = np.array(yg)

    np.save('./Data/synthetic/error_center_stats.npy',{'x0':x0,'y0':y0,'xg':xg,'yg':yg})
    ### answer is the std_dev ~ 0.5 pixels

    bins = 15
    f,axs = plt.subplots(2,2,figsize=(12,10))
    axs = axs.flatten()
    my_hist(axs[0], xg, bins=bins)
    axs[0].set_xlabel('X guess')
    my_hist(axs[1], yg, bins=bins)
    axs[1].set_xlabel('Y guess')
    my_hist(axs[2], x0, bins=bins)
    #axs[2].set_xlim(axs[0].get_xlim())
    axs[2].set_xlabel('X found')
    my_hist(axs[3], y0, bins=bins)
    #axs[3].set_xlim(axs[1].get_xlim())
    axs[3].set_xlabel('Y found')
    f.tight_layout()
    f.savefig('./Data/synthetic/error_center.png')
    plt.show()

def ring_sums():
    data = h5_2_dict('./Data/synthetic/error_test_01.h5')
    A = data['2Ddata'] 
    binsize = 0.01
    
    ### build a common binarr ###
    ri = 2000
    imax = int((ri ** 2. - 2. * ri - 1.)/(1. + 2. * ri) / binsize)
    norm_radius = np.sqrt(2*binsize*ri + binsize**2)
    binarr = np.sqrt(range(1, imax))*norm_radius
    r = np.concatenate(([0.0], binarr))
    r = 0.5 * (r[0:-1] + r[1:])
    ### do the rest of the ringsum stuff ###
    def ringsum2(data, binarr, x0, y0):
        ny, nx = data.shape
        x = 1.*np.arange(0, nx, 1)
        y = 1.*np.arange(0, ny, 1)
        xx,yy = np.meshgrid(x-x0, y-y0)
        R = np.sqrt(xx**2 + yy**2)
        sig,_ = np.histogram(R, bins=np.concatenate((np.array([0.]),binarr)), weights=data)
        return sig
    
    ringsums = np.zeros((100,binarr.size))
    for i in range(100):
        print '{0} of 100'.format(i+1)
        x0 = np.random.normal(3007.5, 0.5)
        y0 = np.random.normal(2007.5, 0.5)
        ringsums[i,:] = ringsum2(A, binarr, x0, y0)

    np.save('./Data/synthetic/error_ringsums_b1e-2.npy',{'ringsums':ringsums, 'binarr':r, 'binsize':binsize})

def plot_ringsums_error(fname='error_ringsums.npy', block=True):
    a = np.load('./Data/synthetic/{0}'.format(fname)).item()
    ringsums = a['ringsums']
    r = a['binarr']
    binsize = a['binsize']

    std_dev = np.std(ringsums, axis=0)
    mean = np.mean(ringsums, axis=0)

    f,axs = plt.subplots(2,1,sharex=True,figsize=(10,7))
    axs[0].fill_between(r, mean-std_dev, mean+std_dev, color='r', alpha=0.6, label=r'1$\sigma$')
    axs[0].fill_between(r, mean-2.*std_dev, mean-std_dev, color='g', alpha=0.6, label=r'2$\sigma$')
    axs[0].fill_between(r, mean+std_dev, mean+2.*std_dev, color='g', alpha=0.6)
    axs[0].fill_between(r, mean-3.*std_dev, mean-2.*std_dev, color='b', alpha=0.6, label=r'$3\sigma$')
    axs[0].fill_between(r, mean+2.*std_dev, mean+3.*std_dev, color='b', alpha=0.6)
    axs[0].plot(r, mean, color='k',label='mean') 
    axs[0].set_title('Ringsum error (N=1000) \n from mis-centering with a binsize={0}'.format(binsize),fontsize=22)
    axs[1].set_xlabel('R (px)',fontsize=16)
    axs[0].legend(fontsize=18)
    axs[1].plot(r, 100.*(std_dev/mean))
    ix = np.abs(r-1000).argmin()
    per_err = np.mean(std_dev[ix::]/mean[ix::])
    axs[1].axhline(per_err*100,color='k',lw=2,linestyle='--')
    axs[1].set_title('Percent Error', fontsize=18)
    axs[1].set_ylabel('%',fontsize=16)
    for axx in axs:
        for loc in ['top','right','bottom','left']:
            axx.spines[loc].set_linewidth(2)
        for tick in axx.xaxis.get_major_ticks():
            tick.label1.set_fontsize(16)
        for tick in axx.yaxis.get_major_ticks():
            tick.label1.set_fontsize(16)
        axx.xaxis.set_tick_params(width=3, length=8)
        axx.yaxis.set_tick_params(width=3, length=8)
    plt.show(block=block)

if __name__ == "__main__":
    ring_sums()
    plot_ringsums_error(block=False)
    plot_ringsums_error(fname='error_ringsums_b1.npy', block=False)
    plot_ringsums_error(fname='error_ringsums_b1e-2.npy')
