import numpy as np
import matplotlib.pyplot as plt
from os.path import join
import argparse
import sys
sys.path.append('../')
from tools.images import get_data
from tools.file_io import h5_2_dict
from core.fitting import find_maximum

def main(folder, big_range=[150.,1000.], atol=1.e-3, sm_range=[600.,700.], saveit=False):
    a = h5_2_dict(join(folder, 'ringsum.h5'))
    data = get_data(a['fname'], color=a['color'])
    
    R,Z = np.meshgrid(np.arange(data.shape[1],dtype='float64'),np.arange(data.shape[0],dtype='float64'))
    R -= a['center'][0]
    Z -= a['center'][1]
    
    Theta = -1.*(np.arctan2(Z,-R)-np.pi)
    Rho = np.sqrt(R**2+Z**2)
    del R,Z
    
    f = plt.figure(figsize=(16,9))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2)
    ax3 = plt.subplot2grid((2,2), (0,1))
    
    ax1.imshow(data, cmap='Greys_r', origin='lower')
    ax1.set_aspect('equal')
    ax1.set_xlim([0, data.shape[1]-1])
    ax1.set_ylim([0, data.shape[0]-1])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title('Ringsum Comparison for \n {0}'.format(folder),fontsize=18)
    
    Rho = Rho.flatten()
    Theta = Theta.flatten()
    data = data.flatten()
    big_ix = np.where(np.logical_and(Rho < big_range[1],Rho > big_range[0]))
    data = data[big_ix]
    Rho = Rho[big_ix]
    Theta = Theta[big_ix]
    
    thetas = np.arange(0.,360.)
    sm_thetas = np.arange(0., 360., 45.)
    for i,t in enumerate(sm_thetas):
        ixs = np.isclose(Theta, np.ones_like(Theta)*t*(np.pi/180.), atol=atol)
        rr = Rho[ixs]
        dd = data[ixs]
        ix = np.argsort(rr)
        ax2.plot(rr[ix],dd[ix],label='{0}'.format(t))
        ax1.plot([a['center'][0],a['center'][0]+2000.*np.cos(t*(np.pi/180.))],[a['center'][1],a['center'][1]+2000.*np.sin(t*(np.pi/180.))],lw=2)
    
    big_ix = np.where(np.logical_and(Rho < sm_range[1],Rho > sm_range[0]))
    data = data[big_ix]
    Rho = Rho[big_ix]
    Theta = Theta[big_ix]
    
    pk = np.zeros_like(thetas)
    for i,t in enumerate(thetas):
        ixs = np.isclose(Theta, np.ones_like(Theta)*t*(np.pi/180.), atol=atol)
        rr = Rho[ixs]
        dd = data[ixs]
        ix = np.argsort(rr)
        pk[i] = find_maximum(rr[ix],dd[ix])
   
    ax2.legend(fontsize=14,frameon=False)
    ax2.set_xlim(big_range)
    ax2.set_xlabel('R (px)',fontsize=14)
    
    ix1 = np.abs(a['r']-sm_range[0]).argmin()
    ix2 = np.abs(a['r']-sm_range[1]).argmin()
    rpk = find_maximum(a['r'][ix1:ix2],a['sig'][ix1:ix2])
   
    dev = np.sqrt(np.nanmean(np.abs(pk - rpk)**2))
    

    ax3.plot(thetas,pk, label='peaks from different angles, std_dev={0:.3f} px'.format(dev))
    ax3.axhline(rpk, label='ringsummed peak',color='r',lw=2)
    ax3.legend(fontsize=14,frameon=False)
    ax3.set_xlabel(r'$\theta$ (deg.)',fontsize=14)
    ax3.set_ylabel('Pk location (pixels)',fontsize=14)
    ax3.set_ylim(sm_range)
    
    for ax in [ax2, ax3]:
        for loc in ['top','right','bottom','left']:
                ax.spines[loc].set_linewidth(1.5)
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(14)
        ax.xaxis.set_tick_params(width=2, length=6)
        ax.yaxis.set_tick_params(width=2, length=6)
    
    plt.tight_layout()
    if saveit:
        plt.savefig(join(folder,'plots/','ringsum_comparison.png'),bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='checks the ringsum regularity for Lyon data')
    parser.add_argument('folder', type=str, help='folder of ringsum.h5 file')
    parser.add_argument('--big_range','-b',type=float, nargs=2, default=[150.,1000.], help='Raidus range to plot')
    parser.add_argument('--pk_range','-p',type=float, nargs=2, default=[600.,700.], help='Radius range to find peak to compare angles')
    parser.add_argument('--atol', type=float, default=1.e-3, help='Absolute tolerance (in radians) of isclose for angles')
    parser.add_argument('--saveit','-s',action='store_true', help='Flag to save figure in plots as ringsum_compare.png')
    args = parser.parse_args()
    
    main(args.folder, big_range=args.big_range, sm_range=args.pk_range, atol=args.atol, saveit=args.saveit)
