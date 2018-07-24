import numpy as np
import matplotlib.pyplot as plt

def lande_g(L,S,J):
    '''
    calculates the lande g factor given L, S and J
    
    g = 1 + {J(J+1) + S(S+1) - L(L+1)} / {2J(J+1)}
    
    Reference: Sobel'man, I.I. Introduction to the 
        Theory of Atomic Spectra. 1972. pp. 277

    Args:
        L (float): L number
        S (float): S number
        J (float): J number
    Returns:
        g (float): lande g factor
    '''
    return 1. + (J*(J+1.) + S*(S+1.) - L*(L+1.)) / (2.*J*(J+1.))

def zeeman_factors(L_up, S_up, J_up, L_lo, S_lo, J_lo, parallel=True):
    '''
    calculates the Zeeman factors and relative amplitudes for a given
    line transition from (L,S,J)_upper to (L,S,J)_lower
    
    assumes LS-coupling (not jj) therefore Z<=30 and that this will
    be used for a weak-field (B < 10T) Zeeman calculation

    Reference: Sobel'man, I.I. Introduction the Theory of
        Atomic Spectra. 1972. pp. 277-280
    
    Notes on getting L,S,J for a given spectral line:
    NIST spectral database will return the upper and lower level
    electron configuration information for a given line. The term
    will be of the form:
            (#){S,P,D,F}
    where #==(2S+1), and {S,P,D,F} correspond to L=0,1,2,3
    Also, NIST will return the value of J for the upper and lower level

    Args:
        L_up (float): the upper state L number
        S_up (float): the upper state S number
        J_up (float): the upper state J number
        L_lo (float): the lower state L number
        S_lo (float): the lower state S number
        J_lo (float): the lower state J number
        parallel (bool, default=True): a boolean flag for the direction
            of observation wrt the magnetic field. If true, the spectrum
            is observed on a line of sight parallel with the magnetic
            field and there are no pi-components of the Zeeman split. If
            false the observation is perpendicular to the field and there
            will be pi-components of the Zeeman split.

    Returns:
        out (dict): a dictionary with the following items
            sp_fac (list) -> sigma plus components zeeman factors
            sp_amp (list) -> relative amplitudes of sigma plus components
            sm_fac (list) -> sigma minus components zeeman factors
            sm_amp (list) -> relative amplitudes of sigma minus components
            if parallel is False:
            pi_fac (list) -> pi components zeeman factors
            pi_amp (list) -> relative amplitudes of pi components
    '''
    g_up = lande_g(L_up,S_up,J_up)
    g_lo = lande_g(L_lo,S_lo,J_lo)
    M_up = np.arange(-J_up,J_up+0.1,1.)
    M_lo = np.arange(-J_lo,J_lo+0.1,1.)

    ## Pi when M_up == M_lo ##
    pi_fac = []
    pi_amp = []
    ## Sigma plus when M_up == M_lo - 1 ##
    sp_fac = []
    sp_amp = []
    ## Sigma minus when M_up == M_lo + 1 ##
    sm_fac = []
    sm_amp = []
    for i,m_up in enumerate(M_up):
        for j,m_lo in enumerate(M_lo):
            if m_up == m_lo:
                ### Pi component ###
                if parallel:
                    pi_fac.append(0)
                    pi_amp.append(0)
                else:
                    pi_fac.append(g_up * m_up - g_lo * m_lo)
                    if J_up == J_lo:
                        pi_amp.append(m_up**2)
                    elif J_up == J_lo + 1.:
                        pi_amp.append(J_up**2 - m_up**2)
                    elif J_up == J_lo - 1.:
                        pi_amp.append((J_up+1)**2 - m_up**2)
            elif m_up == m_lo + 1.:
                ### sigma minus component ###
                sm_fac.append(g_up * m_up - g_lo * m_lo)
                if J_up == J_lo:
                    if parallel:
                        sm_amp.append(0.5*(J_up + m_up)*(J_up + 1. - m_up))
                    else:
                        sm_amp.append(0.25*(J_up + m_up)*(J_up + 1. - m_up))
                elif J_up == J_lo + 1.:
                    if parallel:
                        sm_amp.append(0.5*(J_up + m_up)*(J_up - 1. + m_up))
                    else:
                        sm_amp.append(0.25*(J_up + m_up)*(J_up - 1. + m_up))
                elif J_up == J_lo - 1.:
                    if parallel:
                        sm_amp.append(0.5*(J_up + 1. - m_up)*(J_up - m_up + 2.))
                    else:
                        sm_amp.append(0.25*(J_up + 1. - m_up)*(J_up - m_up + 2.))
            elif m_up == m_lo - 1.:
                ### sigma plus component ###
                sp_fac.append(g_up * m_up - g_lo * m_lo)
                if J_up == J_lo:
                    if parallel:
                        sp_amp.append(0.5*(J_up - m_up)*(J_up + 1. + m_up))
                    else:
                        sp_amp.append(0.25*(J_up - m_up)*(J_up + 1. + m_up))
                elif J_up == J_lo + 1.:
                    if parallel:
                        sp_amp.append(0.5*(J_up - m_up)*(J_up - 1. - m_up))
                    else:
                        sp_amp.append(0.25*(J_up - m_up)*(J_up - 1. - m_up))
                elif J_up == J_lo - 1.:
                    if parallel:
                        sp_amp.append(0.5*(J_up + 1. + m_up)*(J_up + m_up + 2.))
                    else:
                        sp_amp.append(0.25*(J_up + 1. + m_up)*(J_up + m_up + 2.))
    if parallel:
        return {'sp_fac':sp_fac,'sp_amp':sp_amp,'sm_fac':sm_fac,'sm_amp':sm_amp}
    else:
        return {'pi_fac':pi_fac, 'pi_amp':pi_amp, 'sp_fac':sp_fac, 'sp_amp':sp_amp, 'sm_fac':sm_fac, 'sm_amp':sm_amp}

def zeeman_lambda(w0, B, factors, amps=None):
    '''
    calculates the zeeman split shifts from a central wavelength for
    a given magnetic field strength

    Args:
        w0 (float): central wavelength in nm
        B (float): magnetic field strength in T
        factors (list): list of zeeman factors (from zeeman_factors dictionary)
        amps (list, default=None): optional list of amplitudes to renormalize
    Returns:
        lambdas (list): list of zeeman splitting peak locations based on factors
        if amps is not None:
        amps (list): list of normalized amplitudes
    '''
    beta = 9.274e-24 #Bohr magnetron J/Tesla
    lam_hc = 5.034117e15 * w0**2 #lambda^2/hc in nm/J
    lambdas = [x*beta*lam_hc*B + w0 for x in factors]
    if amps is None:
        return lambdas
    else:
        return lambdas, [x/sum(amps) for x in amps]

if __name__ == "__main__":
    print(zeeman_factors(2,0.5,2.5,1,0.5,1.5))
