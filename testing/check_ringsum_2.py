import numpy as np
import matplotlib.pyplot as plt
from os.path import join
import argparse
import sys
sys.path.append('../')
from tools.images import get_data
from tools.file_io import h5_2_dict
from tools.plotting import ringsum_click,ring_plot

def main(folder):
    a = h5_2_dict(join(folder, 'ringsum.h5'))
    data = get_data(a['fname'], color=a['color'])

    R,Z = np.meshgrid(np.arange(data.shape[1],dtype='float64'),np.arange(data.shape[0],dtype='float64'))
    center = a['center']
    R -= center[0]
    Z -= center[1]

    Theta = -1.*(np.arctan2(Z,-R)-np.pi)
    Rho = np.sqrt(R**2+Z**2)
    del R,Z

    Rho = Rho.flatten()
    Theta = np.rad2deg(Theta.flatten())
    Data = data.flatten()

    points,_ = ringsum_click(a['r'],a['sig'],title='pick points')
    binsizes = np.diff(np.concatenate(([0],a['r'])))
    ixs = [np.abs(p - a['r']).argmin() for p in points]
    binsizes = binsizes[ixs]

    ths = []
    dds = []
    for i,p in enumerate(points):
        ixs = np.isclose(Rho,p,atol=binsizes[i])
        datas = Data[ixs]
        thetas = Theta[ixs]
        ixx = np.argsort(thetas)
        dds.append(datas[ixx])
        ths.append(thetas[ixx])

    fig = plt.figure(figsize=(16,9))

    lp = len(points)
    hp = int(lp/2)

    axpic = plt.subplot2grid((lp,2),(0,1),rowspan=lp-hp)
    ring_plot(data,fax=(fig,axpic))
    axpic.plot([center[0],center[0]+2000.],[center[1]]*2,color='white',linestyle='--',lw=2)
    axpic.text(center[0]+1200,center[1]+100,r'$\theta=0$',fontsize=18,color='white')
    axpic.set_xticks([])
    axpic.set_yticks([])

    axbig = plt.subplot2grid((lp,2),(lp-hp,1),rowspan=hp)
    axbig.plot(a['r'],a['sig'],color='k')
    axbig.set_xlim([0,1000])
    axbig.set_xlabel('R (pixels)')
    ylim = [0,12000]
    axbig.set_ylim(ylim)

    axs = []
    for i,p in enumerate(points):
        axs.append(plt.subplot2grid((lp,2),(i,0)))

    for i,p in enumerate(points):
        axs[i].plot(ths[i],dds[i],color='C{0}'.format(i))
        axs[i].set_ylim(ylim)
        axs[i].set_xticklabels([])
        axs[i].set_xticks([0,45,90,135,180,225,270,315])
        axbig.axvline(p,color='C{0}'.format(i),linewidth=2)
    axs[-1].set_xlabel(r'$\theta$ (deg.)')
    axs[-1].set_xticklabels([str(x) for x in [0,45,90,135,180,225,270,315]])
    plt.tight_layout()
    plt.savefig('/home/flanagan/ring_variation.png',bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main('../Data/2017_12_11/0107/')
