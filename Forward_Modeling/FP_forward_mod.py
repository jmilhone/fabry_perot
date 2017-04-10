import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class InputSpec(object):
    def __init__(self, temp=1., vel=0., lam0=488.):
        self.temp = temp    #eV
        self.vel = vel  #km/s
        self.lam0 = lam0    #nm
        if lam0 == 488.:
            self.mu = 40.
        else:
            self.mu = 4.
        self.calc_coeffs()

    def calc_coeffs(self):
        self.wavelength = self.lam0 * (1. - 3.336e-6 * self.vel)  # nm, +vel=blueshift and -vel=redshift
        self.sig = 3.276569e-5 * np.sqrt((self.temp * self.lam0**2)/self.mu)    #nm
        self.FWHM = 2. * np.sqrt(2. * np.log(2.)) * self.sig    #nm
        self.bounds = (self.wavelength - 5.*self.FWHM, self.wavelength + 5.*self.FWHM)

    def plot(self):
        lam = np.linspace(self.bounds[0], self.bounds[1], 1000)
        spec = (2. * np.pi)**(-2) * (self.sig)**(-1) * np.exp(-0.5 * ((lam - self.wavelength) / self.sig)**2)
        f, ax = plt.subplots()
        ax.plot(lam, spec, lw=2)
        ax.plot([self.lam0]*2, [0, spec.max()], 'k--')
        ax.legend(loc='best')
        plt.show()

    def eval_spec(self, lam):
        spec = (2. * np.pi) ** (-2) * (self.sig) ** (-1) * np.exp(-0.5 * ((lam - self.wavelength) / self.sig) ** 2)
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

class OutputFP(object):
    def __init__(self, d, wavelength, f, r):
        self.r = r
        self.costh = f / np.sqrt(f**2 + self.r**2)
        self.m_max = 2. * d/wavelength * 1.e6
        self.m = self.m_max * self.costh
        self.num_pks = int(np.floor(np.floor(self.m_max) - self.m.min() + 1))
        self.pk_r = f * np.sqrt((self.m_max / (np.floor(self.m_max) - np.arange(self.num_pks)))**2 - 1.)

class FP(object):
    def __init__(self, d=0.88, F=20.):
        self.d = d  #mm
        self.F = F  #finesse
        self.Q = (2. / np.pi) * self.F ** 2  # quality factor

    def eval_airy(self, lam, n, costh):
        m = ((2. * self.d * 1.e6)/lam)*costh - n * np.floor(((2. * self.d * 1.e6)/lam))
        return (1. + self.Q * np.sin(np.pi * m)**2)**(-1)

    def convolve_quad(self, npks, costh):
        out = np.zeros((npks, costh.size))
        for i,n in enumerate(range(npks)):
            for j, c in enumerate(costh):
                print i, j
                out[i, j] = quad(lambda x: self.InputSpec.eval_spec(x) * self.eval_airy(x, n, c), self.InputSpec.bounds[0], self.InputSpec.bounds[1])[0]
        return np.sum(out, axis=0)

    def convolve_trapz(self, costh, ns=None):
        lam = np.linspace(self.InputSpec.bounds[0], self.InputSpec.bounds[1], 1000, endpoint=True)
        ll, cth = np.meshgrid(lam, costh, indexing='ij')
        out = np.zeros_like(costh)
        if ns is None:
            ns = self.Output.num_pks
        for i in range(ns+1):
            out += np.trapz(self.InputSpec.eval_spec(ll) * self.eval_airy(ll, i, cth), lam, axis=0)
            print "{0} of {1} done".format(i, ns)
        return out/out.max()

    def run_spec(self, input_spec=InputSpec(), ccd=CCD(), ns=None):
        self.InputSpec = input_spec
        self.CCD = ccd
        self.Output = OutputFP(self.d, self.InputSpec.wavelength, self.CCD.f, self.CCD.r_arr)
        self.instrument_func = self.eval_airy(self.InputSpec.wavelength, 0., self.Output.costh)
        self.Output.output = self.convolve_trapz(self.Output.costh, ns=ns)

    def plot(self):
        f, ax = plt.subplots()
        ax.plot(self.Output.r, self.Output.output, label='Output')
        ax.plot(self.Output.r, self.instrument_func, label='FP Func.')
        for rr in self.Output.pk_r:
            ax.plot([rr]*2, [0, 1], 'k--')
        ax.legend(loc='best')
        plt.show(block=False)

if __name__ == "__main__":
    a = FP()
    a.run_spec()
    a.plot()
    plt.show()