import numpy as np
import matplotlib.pyplot as plt
from os.path import join
import argparse
import sys
sys.path.append('../')
from tools.images import get_data
from tools.file_io import h5_2_dict
from core.fitting import find_maximum,determine_fit_range
from tools.plotting import ringsum_click

def remove_prof(r, sig, max_r=None, poly_num=5):
    '''
    removes a polyfit of order minima from ringsum

    Args:
        r (np.ndarray): binarr of ringsum
        sig (np.ndarray): ringsum to be subtracted
        max_r (list, default=None): r location of peaks, if None
            an interactive click plot will show up
        poly_num (int, default=5): order used for polyfit
    
    Returns:
        max_r (list): peak locations in r
        poff (np.ndarry): fit to be subtracted from sig
    '''
    if max_r is None:
        max_r, _ = ringsum_click(r,sig,title='Please click on peaks')
   
    peak_list = []
    val_list = []
    for pk in max_r:
        idxs = determine_fit_range(r, sig, pk, thres=0.1)
        pk_r,pk_val = find_maximum(r[idxs],sig[idxs],returnval=True)
        peak_list.append(pk_r)
        val_list.append(pk_val)

    peak_loc = np.array(peak_list)
    peak_val = np.array(val_list)

    p = np.polyfit(peak_loc,peak_val,poly_num)
    poff = np.polyval(p,r)
    return peak_loc, poff

def main(folder,err=0.1,center=None):
    a = h5_2_dict(join(folder, 'ringsum.h5'))
    data = get_data(a['fname'],a['color'])

    R,Z = np.meshgrid(np.arange(data.shape[1],dtype='float64'),np.arange(data.shape[0],dtype='float64'))
    if center is None:
        center = a['center']
    R -= center[0]
    Z -= center[1]

    Theta = -1.*(np.arctan2(Z,-R)-np.pi)
    Rho = np.sqrt(R**2+Z**2)
    del R,Z

    Rho = Rho.flatten()
    Theta = Theta.flatten()
    data = data.flatten()
    big_range = [a['r'][0],a['r'][-1]]
    big_ix = np.where(np.logical_and(Rho < big_range[1],Rho > big_range[0]))
    Rho = Rho[big_ix]
    Theta = Theta[big_ix]
    data = data[big_ix]
    
    Rs = []
    Ds = []
    thetas = np.arange(0,360.,45)

    for i,t in enumerate(thetas):
        print t
        ixs = np.isclose(Theta, np.ones_like(Theta)*t*(np.pi/180.),atol=1.e-3)
        rr = Rho[ixs]
        dd = data[ixs]
        ix = np.argsort(rr)
        dd = dd[ix]
        rr = rr[ix]
        #_,p = remove_prof(rr,dd)
        dd /= dd.max()
        Rs.append(rr)
        Ds.append(dd)
    
    ringsum_err = err*a['sig'] + err*np.ones_like(a['sig'])*(a['sig'].max()-a['sig'].min())

    f,axs = plt.subplots(figsize=(16,9))
    axs.fill_between(a['r'],a['sig']+ringsum_err,a['sig']-ringsum_err,alpha=0.3)
    for i,r in enumerate(Rs):
        axs.plot(r,Ds[i],label='{0}'.format(thetas[i]))
    axs.set_xlim([450,750])
    axs.legend()
    plt.show()

if __name__ == "__main__":
    main('../Data/2017_12_11/0109/',err=0.15)



