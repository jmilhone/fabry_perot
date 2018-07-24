from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import multiprocessing as mp
import h5py

def airy_func(wavelength, cos_th, d, F):
    Q = (2. * F / np.pi)**2
    airy = 1.0 / (1.0 + Q * np.sin(np.pi * 2.e6 * d * cos_th / wavelength)**2)
    return airy

def doppler_calc(w0, mu, temp, v):
    sigma = w0 * 3.2765e-5 * np.sqrt(temp / mu)
    w = w0 * (1.0 - 3.336e-9 * v)
    return sigma, w

def gaussian(wavelength, w, sigma, amp=1.):
    norm = 1. / (sigma * np.sqrt(2.*np.pi))
    exp = np.exp(-0.5 * (wavelength-w)**2 / sigma**2)
    return amp * norm * exp

def lorentzian(wavelength, w, gamma, amp=1.):
    A = (amp * 0.5 * gamma) / np.pi
    return A / ((wavelength - w)**2 + (0.5 * gamma)**2)

def forward_model(r, L, d, F, w0, mu, amp, temp, v, nlambda=1024, sm_ang=False, nprocs=6):
    if type(w0) in [list, tuple]:
        if not all([type(x) in [list,tuple] for x in [mu, amp, temp, v]]):
            raise ValueError('need to have a list for all spec params')
        if not all([len(x) == len(w0) for x in [mu, amp, temp, v]]):
            raise ValueError('spec params are not all the same length')

        sigma = []
        w = []
        for i,ww in enumerate(w0):
            s, l = doppler_calc(ww, mu[i], temp[i], v[i])
            sigma.append(s)
            w.append(l)
        
        if nprocs > 1:
            wavelength = np.linspace(min(w) - 10.*max(sigma), max(w) + 10.*max(sigma), nlambda)
        else:
            wavelength = np.linspace(min(w) - 10.*max(sigma), max(w) + 10.*max(sigma), nlambda)[:,np.newaxis]
        
        spec = 0.0
        for idx,ww in enumerate(w):
            spec += gaussian(wavelength, ww, sigma[idx], amp[idx])

    else:
        if not all([type(x) not in [list,tuple] for x in [mu, amp, temp, v]]):
            raise ValueError('need to have a list or not for all spec params')
        sigma, w = doppler_calc(w0, mu, temp, v)
        if nprocs > 1:
            wavelength = np.linspace(w - 10.*sigma, w + 10.*sigma, nlambda)
        else:
            wavelength = np.linspace(w - 10.*sigma, w + 10.*sigma, nlambda)[:,np.newaxis]

        spec = gaussian(wavelength, w, sigma, amp)
    if sm_ang:
        cos_th = 1.0 - 0.5 * (r/L)**2
    else:
        cos_th = L / np.sqrt(L**2 + r**2)
    
    if nprocs > 1:
        def par_func(cos_th, spec, wavelength, d, F, out=None, label=None):
            model = np.zeros_like(cos_th)
            for i, cth in enumerate(cos_th):
                model[i] = trapz(spec*airy_func(wavelength, cth, d, F), wavelength, axis=0)
            if out and label:
                out.put((label, model))
            else:
                return model

        cos_ths = np.array_split(cos_th, nprocs)
        procs = []
        sigs = {}
        out = mp.Queue()
        labels = ['{0}'.format(x) for x in range(nprocs)]
        for k in range(nprocs):
            p = mp.Process(target=par_func, args=(cos_ths[k], spec, wavelength, d, F), 
                    kwargs={'out':out, 'label': labels[k]})
            procs.append(p)
            p.start()

        for i in range(nprocs):
            tup = out.get()
            sigs[tup[0]] = tup[1]

        for p in procs:
            p.join()
        
        model = []
        for k in labels:
            model.append(sigs[k])
        model = np.concatenate(model)
    else:
        model = trapz(spec*airy_func(wavelength, cth, d, F), wavelength, axis=0)
    return model

def ccd(npx=(4016, 6016),px_size=0.004):
    cntr = [(x-1)/2. for x in npx]
    return px_size * np.fromfunction(lambda i, j: np.sqrt((i-cntr[0])**2 + (j-cntr[1])**2), npx, dtype=float)

def ccd_quad(npx=(4016, 6016), px_size=0.004):
    end = (int(npx[0]/2.), int(npx[1]/2.))
    return px_size * np.fromfunction(lambda i,j: np.sqrt((i+0.5)**2+(j+0.5)**2), end, dtype=float)

def recomb_quad(a):
    b = np.vstack((a[::-1,:],a))
    return np.hstack((b[:,::-1],b))

def full_pattern(L,d,F,w0,mu,T,V,A=1.,sm_ang=False,nprocs=6,plotit=False,saveit=None):
    '''
    produces full synthethic ring pattern for Nikon d5200/5300 camera

    Inputs:
        L (float): camera lens focal length (in mm)
        d (float): etalon spacing (in mm)
        F (float): etalon finesse
        w0 (float or list of floats): wavelengths of spectrum (in nm)
        mu (float or list of floats): atomic mass of elements used, same
            order as w0 list (in a.m.u.)
        V (float or list of floats): flow velocity of spectrum (in km/s)
        A (float or list of floats): amplitudes of spectrum, default is 1
        sm_ang (bool, default=F): flag to use the small angle approximation,
            for true synthetic data this should be False
        nprocs (int, default=6): number of processors to use for calc.
        plotit (bool, default=F): flag to plot the resulting rings
        saveit (str, default=None): hdf5 filename for optional saved rings,
            if left to None, the data will not be saved

    Outputs:
        rings (np.ndarray): output of forward forward model
    '''
    if type(w0) is list and type(A) is float:
        A = [A]*len(w0)

    a = ccd_quad()
    rings = forward_model(a.flatten(),L,d,F,w0,mu,A,T,V,sm_ang=sm_ang,nprocs=nprocs)
    rings = rings.reshape(a.shape)
    rings = recomb_quad(rings)

    r = np.arange(0., (np.sqrt(rings.shape[0]**2 + rings.shape[1]**2)/2.) + 0.0005, 0.001)
    ringsum = forward_model(r*0.004, L, d, F, w0, mu, A, T, V, sm_ang=sm_ang) 
    
    if saveit is not None:
        with h5py.File(saveit,'w') as hf:
            hf.create_dataset('2Ddata',data=rings,compression='lzf')
            hf.create_dataset('1Ddata',data=ringsum,compression='lzf')
            hf.create_dataset('1Dr',data=r,compression='lzf')
            hf.create_dataset('L',data=L)
            hf.create_dataset('d',data=d)
            hf.create_dataset('F',data=F)
            hf.create_dataset('amp',data=A)
            hf.create_dataset('w0',data=w0)
            hf.create_dataset('mu',data=mu)
            hf.create_dataset('temp',data=T)
            hf.create_dataset('vel',data=V)
    
    if plotit:
        f,axs = plt.subplots(figsize=(10,7))
        axs.imshow(rings, cmap='Greys_r', interpolation='nearest', vmin=0, origin='lower')
        axs.set_aspect('equal')
        plt.show(block=False)
        f,axs = plt.subplots(figsize=(10,7))
        axs.plot(r,ringsum,lw=2)
        axs.set_xlabel('R (px)')
        plt.show()
    
    return rings

