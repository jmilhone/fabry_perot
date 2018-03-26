import numpy as np
from tools.images import get_data, get_metadata, check_nef
from tools.plotting import center_plot, ringsum_click, ring_plot
from core.ringsum import smAng_ringsum, locate_center
from tools.file_io import dict_2_h5, prep_folder
import matplotlib.pyplot as plt
from os.path import join, abspath
import argparse
import h5py 
from scipy.stats import norm

def get_ringsum(data,x0,y0,binsize=1.0):
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

def remove_prof(r, sig, min_r=None, poly_num=5):
    '''
    removes a polyfit of order minima from ringsum

    Args:
        r (np.ndarray): binarr of ringsum
        sig (np.ndarray): ringsum to be subtracted
        min_r (list, default=None): r location of minima, if None
            an interactive click plot will show up
        poly_num (int, default=5): order used for polyfit
    
    Returns:
        min_r (list): minimum r locations of minima
        poff (np.ndarry): fit to be subtracted from sig
    '''
    if min_r is None:
        min_r, _ = ringsum_click(r,sig,title='Please click order minima')
    
    min_ix = [np.abs(r-x).argmin() for x in min_r]
    min_s = []
    min_r = []
    ixs = []
    for ix in min_ix:
        try:
            a = sig[ix-5:ix+5].argmin()
            min_s.append(sig[ix-5+a])
            min_r.append(r[ix-5+a])
            ixs.append(ix-5+a)
        except ValueError:
            ixs.append(ix)

    p = np.polyfit(min_r,min_s,poly_num)
    poff = np.polyval(p,r)-sig[ixs[-1]]
    return min_r, poff

def main(fname, bgfname=None, color='b', binsize=0.1, xguess=None, 
        yguess=None, block_center=False, click_center=True, find_center=True,
        sub_prof=False, poly_num=5, plotit=False, write=None, folder='./Data'):

    bgdata = None
    if fname[-2:].lower() == "h5":
       with h5py.File(fname, 'r') as f:
           data = f.get("2Ddata").value
           bgdata = np.abs(norm(scale=0.1).rvs(data.shape))
    else:
        fname = check_nef(fname)
        bgfname = check_nef(bgfname)

        print 'loading image...'
        data = get_data(fname, color=color)

    if find_center:
        if click_center:
            xguess,yguess = center_plot(data)

        x0,y0 = locate_center(data, xguess=xguess, yguess=yguess, 
                block_center=block_center, binsize=binsize, plotit=plotit)

        if plotit:
            fig,ax = plt.subplots(figsize=(10,8))
            ring_plot(data, fax=(fig,ax))
            ax.axhline(y0,color='b')
            ax.axvline(x0,color='r')
            plt.show()
    else:
        if xguess is None:
            x0 = data.shape[1]/2.
        else:
            x0 = xguess
        if yguess is None:
            y0 = data.shape[0]/2.
        else:
            y0 = yguess
    # this is from my annulus finder on a real image, ignore
    #print 'overwrite x0, y0...'
    #x0 =  2923.499075195907
    #y0 =  1984.5553261036002

    print 'performing ringsums...'
    smooth_r, smooth_sig0 = get_ringsum(data, x0, y0, binsize=1.0)
    r, sig0 = get_ringsum(data, x0, y0, binsize=binsize)

    if bgfname is not None or fname[-2:].lower() == "h5":
        print 'removing background...'
        if bgdata is None:
            bgdata = get_data(bgfname, color=color)
        _,smooth_bg = get_ringsum(bgdata, x0, y0, binsize=1.0)
        _,bg = get_ringsum(bgdata, x0, y0, binsize=binsize)
        smooth_sig = smooth_sig0 - smooth_bg
        sig = sig0 - bg
    else:
        smooth_sig = smooth_sig0
        sig = sig0

    if sub_prof:
        print 'subtracting profile fit...'
        min_r,smooth_p = remove_prof(smooth_r, smooth_sig, poly_num=poly_num)
        _,p = remove_prof(r, sig, min_r=min_r, poly_num = poly_num)
        smooth_sig -= smooth_p
        sig -= p

    if fname[-2:].lower() == "h5": 
        dic = {'fname': abspath(fname), 'color': color, 'center': (x0, y0),
               'smooth_r': smooth_r, "smooth_sig": smooth_sig, 
               'binsize': binsize, 'r': r, 'sig': sig}
    else:
        a = get_metadata(fname)
        dic = {'date':a['date'], 'time':a['time'], 'fname':abspath(fname),
                    'color':color, 'center':(x0,y0),
                    'smooth_r':smooth_r, 'smooth_sig':smooth_sig, 
                    'binsize':binsize, 'r':r, 'sig':sig}
    if bgfname is not None:
        dic['bg_fname'] = abspath(bgfname)

    print 'done!'

    if write is not None:
        prep_folder(write)
        dict_2_h5(join(write,'ringsum.h5'),dic)
        
    if write is not None or plotit:
        fig = plt.figure(figsize=(12,8))
        ax0 = fig.add_subplot(2,2,1)
        ax1 = fig.add_subplot(2,2,2,sharex=ax0)
        ax2 = fig.add_subplot(2,2,3,sharex=ax0,sharey=ax0)
        ax3 = fig.add_subplot(2,2,4,sharex=ax0,sharey=ax1)

        ax0.plot(smooth_r, smooth_sig, color='r',lw=2,zorder=10)
        ax0.axhline(0,color='k',alpha=0.7)
        ax0.set_title('binsize=1')
        ax1.plot(r, sig, color='r',lw=2,zorder=10)
        ax1.axhline(0,color='k',alpha=0.7)
        ax1.set_title('binsize={0}'.format(binsize))

        ax2.plot(smooth_r, smooth_sig0, lw=2,label='raw',zorder=10)
        ax3.plot(r, sig0, lw=2,label='raw',zorder=10)
        ax2.axhline(0,color='k',alpha=0.7)
        ax3.axhline(0,color='k',alpha=0.7)
        if bgfname is not None:
            ax2.plot(smooth_r, smooth_bg, lw=2,label='background',zorder=10)
            ax3.plot(r, bg, lw=2,label='background',zorder=10)
        if sub_prof:
            ax2.plot(smooth_r, smooth_p, lw=2,label='intensity profile',zorder=10)
            ax3.plot(r, p, lw=2,label='intensity profile',zorder=10)
            for r in min_r:
                ax0.axvline(r,color='k',ls='--',zorder=1,alpha=0.7)
                ax1.axvline(r,color='k',ls='--',zorder=1,alpha=0.7)
                ax2.axvline(r,color='k',ls='--',zorder=1,alpha=0.7)
                ax3.axvline(r,color='k',ls='--',zorder=1,alpha=0.7)
        ax2.set_xlabel('R (px)')
        ax3.set_xlabel('R (px)')
        ax2.legend()
        ax3.legend()
        if write:
            prep_folder(join(write,'plots'))
            plt.savefig(join(write,'plots','ringsum.png'),bbox_inches='tight')
        if plotit:
            plt.show()
        else:
            plt.close()

    return dic

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='stores image ringsum into Data folder')
    parser.add_argument('fname',type=str,help='NEF filename of image you wish to process')
    parser.add_argument('--background','-bg',type=str,default=None,help='NEF filename of\
            background image')
    parser.add_argument('--color', '-c', type=str, default='b', help='r,g,b color to use\
            when ringsumming image: default is b')
    parser.add_argument('--binsize', '-b', type=float, default=0.1, help='binsize for fine\
            ringsum: default is 0.1')
    parser.add_argument('-xy', type=float, nargs='+', default=None, help='x and y guess for\
            the center of image, if not provided the center of the ccd is assumed')
    parser.add_argument('--no_search', action='store_true',help='supresses the center finding\
            algorthim and just uses the provided xy')
    parser.add_argument('--no_click', action='store_true', help='supresses the interactive\
            click plot for getting a center guess')
    parser.add_argument('--block', '-bk', action='store_true', help='sets a 600x600 square\
            block centered on xyguess for center finding')
    parser.add_argument('--sub_prof', '-sp', action='store_true', help='perform a polyfit\
            profile subtraction, launches an interactive plot')
    parser.add_argument('--poly_num', type=int, default=5, help='polyfit order to use for\
            profile subtraction: default is 5')
    parser.add_argument('--no_plot', action='store_true', help='supresses plots')
    parser.add_argument('--write', '-w', type=str, default=None, help='saves data and plot to\
            specified folder, default is None which will not write data.')
    args = parser.parse_args()

    plotit = not args.no_plot
    click_center = not args.no_click
    find_center = not args.no_search

    if args.xy is None:
        xguess = None
        yguess = None
    else:
        xguess, yguess = args.xy

    _ = main(args.fname, bgfname=args.background, color=args.color, binsize=args.binsize,
            block_center=args.block, sub_prof=args.sub_prof, poly_num=args.poly_num,
            write=args.write, plotit=plotit, click_center=click_center, xguess=xguess,
            yguess=yguess, find_center=find_center)

