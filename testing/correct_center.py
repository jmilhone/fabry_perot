import numpy as np
import matplotlib.pyplot as plt
from os.path import join
import argparse
import sys
import time
from scipy.optimize import fmin
import multiprocessing as mp
sys.path.append('../')
from tools.images import get_data
from tools.file_io import h5_2_dict
from core.fitting import find_maximum

def _par_func(thetas, Theta, Rho, data, out=None, label=None):
    pk = np.zeros_like(thetas)
    for i,t in enumerate(thetas):
        ixs = np.isclose(Theta, np.ones_like(Theta)*t*(np.pi/180.), atol=1.e-3)
        rr, dd = Rho[ixs], data[ixs]
        ix = np.argsort(rr)
        pk[i] = find_maximum(rr[ix], dd[ix])
    if out and label:
        out.put((label, pk))
    else:
        return pk

def compute_std_dev((x0,y0),data,R,Z,pk_range):
    Theta = -1. * (np.arctan2(Z-y0, -(R-x0))-np.pi)
    Rho = np.sqrt((R-x0)**2+(Z-y0)**2)
   
    big_ix = np.where(np.logical_and(Rho > pk_range[0],Rho < pk_range[1]))
    data = data[big_ix]
    Rho = Rho[big_ix]
    Theta = Theta[big_ix]

    nprocs = 30

    thetas = np.array_split(np.arange(0.,360.),nprocs)
    procs = []
    sigs = {}
    out = mp.Queue()
    labels = ['{0}'.format(x) for x in range(nprocs)]
    for k in range(nprocs):
        p = mp.Process(target=_par_func, args=(thetas[k], Theta, Rho, data), kwargs={'out':out, 'label':labels[k]})
        procs.append(p)
        p.start()

    for i in range(nprocs):
        tup = out.get()
        sigs[tup[0]] = tup[1]

    for p in procs:
        p.join()

    pks = []
    for k in labels:
        pks.append(sigs[k])
    pks = np.concatenate(pks)

    return np.nanstd(pks)

def main(folder, pk_range=[600.,700.]):

    a = h5_2_dict(join(folder, 'ringsum.h5'))
    data = get_data(a['fname'], color=a['color'])

    R,Z = np.meshgrid(np.arange(data.shape[1],dtype='float64'),np.arange(data.shape[0],dtype='float64'))

    x0 = a['center'][0]
    y0 = a['center'][1]

    data = data.flatten()
    R = R.flatten()
    Z = Z.flatten()
    
    print 'starting with x0={0} and y0={1}'.format(x0,y0)
    t0 = time.time()
    xopt = fmin(compute_std_dev, (x0, y0), args=(data, R, Z, pk_range,), disp=True, maxiter=100)
    print 'evaluation completed in {0} seconds'.format(time.time()-t0)
    print xopt


if __name__ == "__main__":
    main('../Data/2017_12_11/0199')
