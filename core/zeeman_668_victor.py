import numpy as np
import matplotlib.pyplot as plt

def lande_g(L,S,J):
    return 1. + (J*(J+1.) + S*(S+1.) - L*(L+1.)) / (2.*J*(J+1.))

def zeeman_factors(L_up,S_up,J_up,L_lo,S_lo,J_lo):
    g_up = lande_g(L_up,S_up,J_up)
    g_lo = lande_g(L_lo,S_lo,J_lo)
    M_up = np.arange(-J_up,J_up+0.1,1.)
    M_lo = np.arange(-J_lo,J_lo+0.1,1.)

    ## Pi when M_up-M_lo==0 ##
    pi_fac = []
    pi_amp = []
    ## Sigma plus when M_up-M_lo==1 ##
    sp_fac = []
    sp_amp = []
    ## Sigma minus when M_up-M_low==-1 ##
    sm_fac = []
    sm_amp = []
    for i,m_up in enumerate(M_up):
        for j,m_lo in enumerate(M_lo):
            if m_up - m_lo == 0:
                ### Pi ###
                pi_fac.append(g_up * m_up - g_lo * m_lo)
                if J_up == J_lo:
                    pi_amp.append(0)
                else:
                    pi_amp.append(0)
            elif m_up - m_lo == 1.:
                ### Sigma Plus ###
                sp_fac.append(g_up * m_up - g_lo * m_lo)
                if J_up == J_lo:
                    sp_amp.append((J_up + m_up)*(J_up + 1. - m_up))
                elif (J_lo - J_up) == -1:
                    sp_amp.append((J_up + m_up)*(J_up - 1. + m_up))
                elif (J_lo - J_up) == +1:
                    sp_amp.append((J_up + 1.- m_up)*(J_up - m_up + 2))
            elif m_up - m_lo == -1.:
                ### Sigma Minus ###
                sm_fac.append(g_up * m_up - g_lo * m_lo)
                if J_up == J_lo:
                    sm_amp.append((J_up - m_up)*(J_up + 1. + m_up))
                elif (J_lo-J_up) == -1:
                    sm_amp.append((J_up - m_up)*(J_up - 1. - m_up))
                elif (J_lo-J_up) == +1:
                    sm_amp.append((J_up +1. + m_up)*(J_up + m_up + 2.))
    
    out = {'pi_fac':pi_fac,'pi_amp':pi_amp, 'sp_fac':sp_fac,
            'sp_amp':sp_amp,'sm_fac':sm_fac,'sm_amp':sm_amp}
    return out

def zeeman_lambda(w0, B, factors, amps=None):
     
    beta = 9.274e-24 #Bohr magnetron J/Tesla
    lam_hc = 5.034117e15 * w0**2 #lambda^2/hc in nm/J
    lambdas = [x*beta*lam_hc*B + w0 for x in factors]
    if amps is None:
        return lambdas
    else:
        return lambdas, [x/sum(amps) for x in amps]

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
        amp (float, default=1.0): amplitude of spectrum, if 1 then the
                spectrum will integrate to 1 over infinity

    Returns:
        spec (np.ndarray): spectrum evaluated on wavelength array
    '''
    if norm:
        norm = 1. / (sigma * np.sqrt(2.*np.pi))
    else:
        norm = 1.
    exp = np.exp(-0.5 * (wavelength-w)**2 / sigma**2)
    return amp * norm * exp

def main(temp,I,plotit=True):
    a = zeeman_factors(2,0.5,2.5,1,0.5,1.5)
    w0 = 487.12
#    a = zeeman_factors(3,1.5,3.5,2,1.5,2.5)
#    w0 = 668.61
    mu = 39.948
    B = (0.0146/80) * I #Tesla
    lambdas, amps = zeeman_lambda(w0, B, a['sm_fac']+a['sp_fac'], amps=a['sm_amp']+a['sp_amp'])
    print amps

    f,ax = plt.subplots(figsize=(10,7))
    ax.axvline(w0,color='k',lw=1,linestyle='--')
    sblah,_ = doppler_calc(w0,mu,temp,0.0)
    mn = w0 - 5*sblah
    mx = w0 + 5*sblah
    wavelength = np.linspace(mn,mx,2048,endpoint=True)
    spec = 0.
    for l,a in zip(lambdas,amps):
        ax.plot([l]*2,[0,a],lw=2,color='lightseagreen')
        sigma, w = doppler_calc(l,mu,temp,0.0)
        spec += gaussian(wavelength, w, sigma, amp=a, norm=False)
        ax.fill_between(wavelength, gaussian(wavelength, w, sigma, amp=a, norm=False), color='b', alpha=0.4)
    ax.fill_between(wavelength, spec, color='r', alpha=0.4)
    ax.set_xlim([mn,mx])
    ax.set_ylim([0,None])
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    if plotit:
        plt.show()
    else:
        plt.close()

    return wavelength, spec

if __name__ == "__main__":
    #wavelength, spec = main(0.8, 600)
    a = zeeman_factors(2,0.5,2.5,1,0.5,1.5)
    print a['sm_fac']+a['sp_fac']
    print a['sm_amp']+a['sp_amp']
    #f,ax = plt.subplots(figsize=(10,7))
    #temps = [0.4,0.6,0.8,1.0,1.2,1.4]
    #for t in temps:
    #    w,s = main(t, 800, plotit=False)
    #    ax.plot(w,s,lw=2,label='I=800A, Ti={0}'.format(t))

    #ax.set_xlabel(r'$\lambda$ (nm)')
    #ax.set_ylabel('amplitude (arb.)')
    #ax.legend(loc='best')
    #ax.get_xaxis().get_major_formatter().set_useOffset(False)
    #ax.set_ylim([0,None])
    #plt.show()
