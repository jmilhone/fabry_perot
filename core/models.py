import numpy as np
from scipy.integrate import trapz
from core.zeeman import zeeman_lambda

def peak_calculator(L, d, w, order):
    '''
    Simple peak calculator for ideal Fabry-Perot.
    
    r_j = L * Sqrt[ ((2d/w)/(floor(2d/w)-j))^2 - 1]

    Args:
        L (float): L value, units will set return units (pixels or mm)
        d (float): etalon spacing in mm
        w (float): peak wavelength in nm
        order (np.ndarray or int): j order number. 0 gives order nearest center
            of rings, increasing j as you move out.

    Returns:
        peaks (np.ndarray or float): the location of the peak(s)
            in the same units as L with the same length as order
    '''
    m = 2.e6 * d / w
    m0 = np.floor(m)
    return L * np.sqrt( m**2 / (m0 - order)**2 - 1.0 )

def airy_func(wavelength, cos_th, d, F):
    '''
    Computes the Airy function (ideal Fabry-Perot instument function)
    as a function of wavelength and cos(theta) with parameters d and F

    Airy = [1. + Q * sin(pi * (2d/w) * cos_th)^2]^(-1)
    Q = (2 * F / pi)^2

    Args:
        wavelength (np.ndarray): wavelength array in nm
        cos_th (np.ndarray): cos(theta) array
        d (float): etalon spacing in mm
        F (float): etalon finesse

    Returns:
        airy (np.ndarray): evaluated airy function
    '''
    Q = (2. * F / np.pi)**2
    airy = 1.0 / (1.0 + Q * np.sin(np.pi * 2.e6 * d * cos_th / wavelength)**2)
    return airy

def doppler_calc(w0, mu, temp, v):
    '''
    Computes the doppler broadening sigma and the new central wavelength
    from the doppler shift
    
    sigma = w0 * Sqrt[kT / mc^2]
    w = w0 * (1 - v/c)
    
    Args: 
        w0 (float): unshifted wavelength in nm
        mu (float): atomic mass in amu
        temp (float): temperature in eV
        v (float): velocity in m/s

    Returns:
        sigma (float): sigma in nm
        w (float): shifted wavelength in nm
    '''
    sigma = w0 * 3.2765e-5 * np.sqrt(temp / mu)
    w = w0 * (1.0 - 3.336e-9 * v)
    return sigma, w

def gaussian(wavelength, w, sigma, amp=1., norm=True):
    '''
    Computes a gaussian for a given central wavelength, sigma and amp

    spec = amp/(sigma*Sqrt[2*pi]) * Exp[ -0.5 * (wavelength - w)^2 / sigma^2 ]
        
    Args:
        wavelength (np.ndarray): wavelength array to calculate spec on
        w (float): central wavelength (same units as wavelength array)
        sigma (float): sigma for gaussian (same units as w)
        amp (float, default=1.0): amplitude of spectrum
        norm (bool, default=True): if true, the gaussian will be normalized
            to integrate to 1 over infinity then the amp factor will be multiplied

    Returns:
        spec (np.ndarray): spectrum evaluated on wavelength array
    '''
    if norm:
        norm = 1. / (sigma * np.sqrt(2.*np.pi))
    else:
        norm = 1.
    exp = np.exp(-0.5 * (wavelength-w)**2 / sigma**2)
    return amp * norm * exp

def lorentzian(wavelength, w, gamma, amp=1.):
    '''
    Computes a lorentzian for a given central wavelength, gamma and amp

    spec = (amp/(2*pi))  * [gamma / ((wavelength-w)**2 + (0.5*gamma)**2)]

    Args:
        wavelength (np.array): wavelength array to calculate spec on
        w (float): central wavelength (same units as wavelength array)
        gamma (float): lorentzian gamma parameter
        amp (float, default=1.): amplitude in addition to one that integrates
            spec to 1 over infinity

    Returns:
        spec (np.ndarray): spectrum evaluated on wavelength array
    '''
    A = (amp * 0.5 * gamma) / np.pi
    return A / ((wavelength - w)**2 + (0.5 * gamma)**2)

def forward_model(r, L, d, F, w0, mu, amp, temp, v, nlambda=1024, sm_ang=True):
    '''
    Convolves the Doppler spectrum with the ideal Fabry-Perot Airy function.

    Args:
        r (np.ndarray): array of r values to compute model on 
        L (float): camera lens focal length, same units as r (pixels or mm)
        d (float): etalon spacing (mm)
        F (float): etalon finesse
        w0 (float or list): central wavelength(s) in nm
        mu (float or list): mass(es) in amu
        temp (float or list): temperature(s) in eV
        v (float or list): velocities in m/s 
        nlambda (int, default=1024): number of points in wavelength array
        sm_ang (bool, default=True): use the small angle approx or not

    Returns:
        model (np.ndarray): array length of r of forward model
    '''
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

        wavelength = np.linspace(min(w) - 10.*max(sigma), 
                max(w) + 10.*max(sigma), nlambda)[:,np.newaxis]
        
        spec = 0.0
        for idx,ww in enumerate(w):
            spec += gaussian(wavelength, ww, sigma[idx], amp[idx])

    else:
        if not all([type(x) not in [list,tuple] for x in [mu, amp, temp, v]]):
            raise ValueError('need to have a list or not for all spec params')

        sigma, w = doppler_calc(w0, mu, temp, v)
        wavelength = np.linspace(w - 10.*sigma, w + 10.*sigma, nlambda)[:,np.newaxis]
        spec = gaussian(wavelength, w, sigma, amp)

    if sm_ang:
        cos_th = 1.0 - 0.5 * (r/L)**2
    else:
        cos_th = L / np.sqrt(L**2 + r**2)

    cos_th = cos_th.reshape((1,len(r)))

    model = trapz(spec*airy_func(wavelength, cos_th, d, F), wavelength, axis=0)
    return model

def match_finesse_forward(r, L, d, F, temp, v, errtemp=None, w0=487.98634, mu=39.948):
    sigma, w = doppler_calc(w0, mu, temp, v*1000.)
    if errtemp is not None:
        errsigma, _ = doppler_calc(w0, mu, errtemp, v*1000.)
        sigma = np.sqrt(sigma**2 + errsigma**2)
    wavelength = np.linspace(w - 10.*sigma, w + 10.*sigma, 512)[:,np.newaxis]
    spec = gaussian(wavelength, w, sigma, norm=False)
    
    cos_th = 1.0 - 0.5 * (r/L)**2
    model = trapz(spec*airy_func(wavelength, cos_th, d, F), wavelength, axis=0)
    return model

def lyon_temp_forward(r, L, d, F, current, T, V, E=None):
    w0 = 487.98634
    mu = 39.948

    ### my previous calculation ###
    #zeeman_fac = [-1.4, -19./15., -17./15., -1., 1., 17./15., 19./15., 1.4]
    #zeeman_amp = [1.0, 3.0, 6.0, 10.0, 6.0, 3.0, 1.0]
    ### Victor's calculation ###
    zeeman_fac = [-1.,-17./15.,-19./15.,-1.4,1.4,19./15.,17./15.,1.]
    zeeman_amp = [20., 12., 6., 2., 2., 6., 12., 20.] 
    
    B = (0.0146/80.) * current
    
    sblah,w = doppler_calc(w0, mu, T, V*1.e3)
    if E is not None:
        eblah,_ = doppler_calc(w0, mu, E, V*1.e3)
        sblah = np.sqrt(sblah**2 + eblah**2)
    lambdas, amps = zeeman_lambda(w, B, zeeman_fac, amps=zeeman_amp)
    mn = w - 10.*sblah
    mx = w + 10.*sblah
    wavelength = np.linspace(mn,mx,1024,endpoint=True)[:,np.newaxis]
    spec = 0.
    for l,a in zip(lambdas,amps):
        sigma, _ = doppler_calc(l, mu, T, 0.0) ## already did the velocity shift
        if E is not None:
            esigma,_ = doppler_calc(l, mu, E, 0.0)
            sigma = np.sqrt(sigma**2 + esigma**2)
        spec += gaussian(wavelength, l, sigma, amp=a, norm=False)

    cos_th = 1.0 - 0.5 * (r/L)**2
    model = trapz(spec*airy_func(wavelength, cos_th, d, F), wavelength, axis=0)
    return model
