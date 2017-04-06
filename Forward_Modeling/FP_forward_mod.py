import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve

class InputSpec(object):
    def __init__(self, temp=1., vel=0., wavelength=488.):
        self.temp = temp    #eV
        self.vel = vel  #km/s
        self.wavelength = wavelength    #nm
        if wavelength == 488.:
            self.mu = 40.
        else:
            self.mu = 4.
        self.calc_coeffs()

    def calc_coeffs(self):
        self.lam0 = self.wavelength * (1. - 3.336e-6 * self.vel)  # nm, +vel=blueshift and -vel=redshift
        self.sig = 3.276569e-5 * np.sqrt((self.temp * self.lam0**2)/self.mu)    #nm
        self.FWHM = 2. * np.sqrt(2. * np.log(2.)) * self.sig    #nm

    def plot(self):
        lam = np.linspace(self.lam0 - 5. * self.FWHM, self.lam0 + 5. * self.FWHM, 1000)
        spec = (2. * np.pi)**(-2) * (self.sig)**(-1) * np.exp(-0.5 * ((lam - self.lam0) / self.sig)**2)
        f, ax = plt.subplots()
        ax.plot(lam, spec, lw=2)
        ax.plot([self.wavelength]*2, [0, spec.max()], 'k--')
        ax.legend(loc='best')
        plt.show()

    def eval_spec(self, lam):
        spec = (2. * np.pi) ** (-2) * (self.sig) ** (-1) * np.exp(-0.5 * ((lam - self.lam0) / self.sig) ** 2)
        return spec

class CCD(object):
    def __init__(self, size=15.6, npx=4096, f=150., binsize=0.001):
        self.size = size
        self.npx = npx
        self.f = f
        self.binsize = binsize
        self.calc_r_arr()
    def calc_r_arr(self):
        self.r_arr = np.arange(0.0, self.size/2., self.binsize)
        self.cosTh = self.f / np.sqrt(self.f**2 + self.r_arr**2)


class FP(InputSpec, CCD):
    def __init__(self, temp=1., vel=0., wavelength=488., size=15.6, npx=4096, f=150., binsize=0.001, d=0.88, F=20.):
        InputSpec.__init__(self, temp=temp, vel=vel, wavelength=wavelength)
        CCD.__init__(self, size=size, npx=npx, f=f, binsize=binsize)
        self.d = d  #mm
        self.F = F  #finesse
        self.Q = (2. / np.pi) * self.F ** 2  # quality factor
        self.calc_m0()
        self.spec = self.eval(self.lam)
        self.airy = (1. + self.Q * np.sin(np.pi * self.m_arr)**2.)**(-1.)

    def calc_m0(self):
        self.m_max = 2. * self.d/self.lam0 * 1.e6   # 2d/lambda, 1.e6 is bc d(mm) and lambda(nm)
        self.m_arr = self.m_max * self.cosTh
        self.lam = 2. * self.d/self.m_arr * 1.e6

    def plot(self):
        f, ax = plt.subplots(2,1)
        ax[0].plot(self.r_arr, self.airy)
        ax[0].set_title('FP Instrument Function')
        ax[1].plot(self.lam, self.spec)
        ax[1].set_title('Input Spec')
        plt.show()

if __name__ == "__main__":
    a = FP()
    a.plot()