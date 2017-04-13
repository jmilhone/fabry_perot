import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class Source(object):
    def __init__(self, temp=1., vel=0., lam0=488., mu=40.):
        self.temp = temp    #eV
        self.vel = vel  #km/s
        self.lam0 = lam0    #nm
        self.mu = mu
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
    def __init__(self, f=37500., npx=(6036, 4020)):
        self.npx = npx
        self.f = f
        self.calc_lin_arr()
        self.calc_pix_arr()

    def calc_pix_arr(self):
        self.center = self.npx[0]/2, self.npx[1]/2
        x = np.linspace(0, self.npx[0]-1, self.npx[0]) - self.center[0]
        y = np.linspace(0, self.npx[1]-1, self.npx[1]) - self.center[1]
        X, Y = np.meshgrid(x, y)
        self.pix_arr = np.sqrt(X**2 + Y**2)
        self.center = (self.center[1], self.center[0])

    def calc_lin_arr(self):
        self.lin_arr = np.linspace(0.0, max(self.npx)/2., 2.*max(self.npx))
        self.cos_th = self.f / np.sqrt(self.f**2 + self.lin_arr**2)

class Etalon(object):
    def __init__(self, d=0.88, F=20.):
        self.d = d  #mm
        self.F = F  #finesse
        self.Q = (2. / np.pi) * self.F ** 2  # quality factor

    def eval_airy(self, lam, cos_th):
        return (1. + self.Q * np.sin(np.pi * 2.e6 * self.d * (1./lam) * cos_th)**2)**(-1)

def eval_lin_fabry(source, ccd, etalon):
    lam_arr = np.linspace(source.bounds[0], source.bounds[1], 1000, endpoint=True)
    ll, cth = np.meshgrid(lam_arr, ccd.cos_th, indexing='ij')
    lin_out = np.trapz(source.eval_spec(ll) * etalon.eval_airy(ll, cth), lam_arr, axis=0)
    return lin_out/lin_out.max()

def interpolate_fabry_to_ccd(lin_arr, lin_data, ccd_grid):
    return griddata(lin_arr, lin_data, ccd_grid, fill_value=0.0)

def main_argon(temp=1.0, vel=0.0, npx=(6036, 4020), f=37500., d=0.88, F=20., plotit=True):
    source = Source(temp=temp, vel=vel, lam0=487.8733)
    ccd = CCD(f=f, npx=npx)
    etalon = Etalon(d=d, F=F)

    inputs = {"T": temp, "V": vel, "npx": npx, "f": f, "d": d, "F": F}

    pixel_arr = ccd.lin_arr
    costh_arr = ccd.cos_th
    print 'Evaluating convolution...'
    lin_data = eval_lin_fabry(source, ccd, etalon)
    instrument_func = etalon.eval_airy(source.wavelength, costh_arr)

    pixel_grid = ccd.pix_arr
    print 'Interpolating to grid...'
    ccd_data = interpolate_fabry_to_ccd(pixel_arr, lin_data, pixel_grid)

    m_max = 2.e6 * d / source.wavelength
    m_arr = m_max * costh_arr
    num_pks = int(np.floor(np.floor(m_max) - m_arr.min() + 1))
    pk_locations = ccd.f * np.sqrt((m_max / (np.floor(m_max) - np.arange(num_pks)))**2 - 1.)
    pk_m_arr = m_max * ccd.f/np.sqrt(ccd.f**2 + pk_locations**2)

    if plotit:
        print 'Plotting...'
        f, ax = plt.subplots()
        ax.plot(pixel_arr, lin_data, label='data')
        ax.plot(pixel_arr, instrument_func, label='instr. f')
        ax.set_xlabel('pixels')
        ax.set_xlim([0, min(npx) / 2.])
        last_pk_ix = np.where(pk_locations > min(npx) / 2.)[0][0]
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(pk_locations[0:last_pk_ix])
        ax2.set_xticklabels(pk_m_arr[0:last_pk_ix].astype(str))
        ax2.set_xlabel(r'm # using $\lambda=${0:.4f}'.format(source.wavelength))
        for loc in pk_locations:
            ax.plot([loc]*2, [0, 1], 'k--')
        ax.legend(loc='best')
        plt.show(block=False)

        f, ax = plt.subplots()
        ax.imshow(ccd_data, cmap='gray', origin='lower')
        ax.set_aspect(1.0)
        plt.show()

    print 'Done!'

    return {"pixel_arr": pixel_arr, "costh_arr": costh_arr, "lin_data": lin_data, "instrument_func": instrument_func,
            "pixel_grid": pixel_grid, "ccd_data": ccd_data, "m_max": m_max, "m_arr": m_arr,
            "pk_locations": pk_locations, "pk_m_arr": pk_m_arr, "inputs": inputs}

if __name__ == "__main__":
    a = main_argon()
