import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import h5py
import cPickle as pickle
import time
import argparse
import model 
class Source(object):
    def __init__(self, temp=1., vel=0., lam0=488., mu=40., amp=None, name='Source generic'):
        self.name = name
        if type(temp) in [float, int]:
            self.temp = temp    #eV
        else:
            self.temp = np.array(temp)
        if type(vel) in [float, int]:
            self.vel = vel  #km/s
        else:
            self.vel = np.array(vel)
        if type(lam0) in [float, int]:
            self.lam0 = lam0    #nm
        else:
            self.lam0 = np.array(lam0)
        if type(mu) in [float, int]:
            self.mu = mu
        else:
            self.mu = np.array(mu)
        if amp is not None:
            if type(amp) in [float, int]:
                self.amp = amp
            else:
                self.amp = np.array(amp)
        else:
            self.amp = amp

        self.calc_coeffs()

    def calc_coeffs(self):
        self.wavelength = self.lam0 * (1. - 3.336e-6 * self.vel)  # nm, +vel=blueshift and -vel=redshift
        self.sig = 3.276569e-5 * np.sqrt((self.temp * self.lam0**2)/self.mu)    #nm
        self.FWHM = 2. * np.sqrt(2. * np.log(2.)) * self.sig    #nm
        if self.amp is None:
            self.amp = (2. * np.pi) ** (-0.5) * (self.sig) ** (-1)
            self.A = None
        else:
            self.A = self.amp
        if type(self.wavelength) is float:
            self.bounds = (self.wavelength - 5.*self.FWHM, self.wavelength + 5.*self.FWHM)
        else:
            ix = np.argsort(self.wavelength)
            wv = self.wavelength[ix]
            fwhm = self.FWHM[ix]
            self.bounds = (wv[0] - 5.*fwhm[0], wv[-1] + 5.*fwhm[-1])

    def plot(self):
        lam = np.linspace(self.bounds[0], self.bounds[1], 1000)
        spec = self.eval_spec(lam)
        f, ax = plt.subplots()
        ax.plot(lam, spec, lw=2)
        if type(self.lam0) is float:
            ax.plot([self.lam0]*2, [0, spec.max()], 'k--')
        plt.show()

    def eval_spec(self, lam):
        if type(self.lam0) is float:
            return self.amp * np.exp(-0.5 * ((lam - self.wavelength) / self.sig) ** 2)
        else:
            out = np.zeros_like(np.array(lam))
            for i, w in enumerate(self.wavelength):
                out += self.amp[i] * np.exp(-0.5 * ((lam - w) / self.sig[i]) ** 2)
            return out

    def params(self):
        return {"lam_0": self.lam0, "temp": self.temp, "vel": self.vel, "mu": self.mu, "A": self.A, "type": self.name}

class Argon(Source):
    def __init__(self, temp=1., vel=0.):
        super(Argon, self).__init__(temp=temp, vel=vel, lam0=487.98634, mu=39.948, name='Argon')

class Mercury(Source):
    def __init__(self, temp=1., vel=0.):
        #super(Mercury, self).__init__(temp=temp, vel=vel, lam0=546.07498, mu=200.59, name='Argon')
        super(Mercury, self).__init__(temp=temp, vel=vel, lam0=576.96095, mu=200.59, name='Argon')

class Thorium(Source):
    def __init__(self, argon_temp=None, argon_amp=None):
        if argon_temp is None and argon_amp is None:
            #temp=1000K -> 0.0861733eV
            super(Thorium, self).__init__(temp=0.0861733, vel=0.0, lam0=487.8733020, mu=232.03806, name='Thorium')
        else:
            if argon_temp is None:
                ar_t = 1.
            else:
                ar_t = argon_temp
            if argon_amp is None:
                ar_amp = 1.
            else:
                ar_amp = argon_amp
            super(Thorium, self).__init__(temp=[0.0861733, ar_t], vel=0.0, lam0=[487.8733020, 487.98634],
                                          mu=[232.03806, 39.948], amp=[1, ar_amp], name='Thorium+Argon')

class Helium(Source):
    def __init__(self, temp=0.5, vel=0., rel_amps=None):
        self.temp = temp
        self.vel = vel
        self.mu = 4.002602
        self.lam0 = np.array([468.5376849, 468.54072254, 468.55244041, 468.55680062, 468.57038494, 468.5704380,
                              468.575707958, 468.5757974, 468.58040922, 468.58308897, 468.588412282, 468.59055531,
                              468.591788438])
        if rel_amps is None:
            #From Evan's Thesis (see Cooper's he4686tivoigtfit.pro)
            rel_amps = [11., 9., 13., 3., 6., 16., 1., 2., 9., 1.5, 1., 28., 1.]
        rel_amps = np.array(rel_amps)
        self.amp = rel_amps/rel_amps.max()
        super(Helium, self).__init__(temp=self.temp, vel=self.vel, lam0=self.lam0, mu=self.mu,
                                     amp=self.amp, name='Helium')

class CCD(object):
    def __init__(self, lens=150., px_size=4., npx=(6036, 4020), lin_pts=5e4, name='CCD generic'):
        self.name = name
        self.npx = npx
        self.lens = lens
        self.px_size = px_size
        self.f =  lens / (px_size * 1.e-3)
        self.calc_lin_arr(lin_pts)
        self.calc_pix_arr()

    def calc_pix_arr(self):
        self.center = self.npx[0]/2, self.npx[1]/2
        x = np.linspace(0, self.npx[0]-1, self.npx[0]) - self.center[0]
        y = np.linspace(0, self.npx[1]-1, self.npx[1]) - self.center[1]
        X, Y = np.meshgrid(x, y)
        self.pix_arr = np.sqrt(X**2 + Y**2)
        self.center = (self.center[1], self.center[0])

    def calc_lin_arr(self, lin_pts):
        self.lin_arr = np.linspace(0.0, max(self.npx)/2., lin_pts, endpoint=True)
        self.cos_th = self.f / np.sqrt(self.f**2 + self.lin_arr**2)

    def params(self):
        return {"type": self.name, "f": (self.lens, self.f), "px_size": self.px_size,
                "npx": self.npx, "center": self.center}

    def arrays(self):
        return {"r": self.lin_arr, "cos_th": self.cos_th, "grid": self.pix_arr}

class NikonD5200(CCD):
    def __init__(self, lens=150., lin_pts=5e4):
        super(NikonD5200, self).__init__(lens=lens, px_size=4.0, npx=(6036, 4020), lin_pts=lin_pts, name='Nikon D5200')

class Andor(CCD):
    def __init__(self, lens=150., lin_pts=5e4):
        super(Andor, self).__init__(lens=lens, px_size=13.0, npx=(1024, 1024), lin_pts=lin_pts, name='Andor iStar334T')

class Etalon(object):
    def __init__(self, d=0.88, F=20.):
        self.d = d  #mm
        self.F = F  #finesse
        self.Q = (2. / np.pi) * self.F ** 2  # quality factor

    def eval_airy(self, lam, cos_th):
        return (1. + self.Q * np.sin(np.pi * 2.e6 * self.d * (1./lam) * cos_th)**2)**(-1)

    def params(self):
        return {"d": self.d, "F": self.F, "Q": self.Q}

# def eval_lin_fabry(source, ccd, etalon):
#     lam_arr = np.linspace(source.bounds[0], source.bounds[1], 1000, endpoint=True)
#     if ccd.cos_th.size < 1.e5:
#         ll, cth = np.meshgrid(lam_arr, ccd.cos_th, indexing='ij')
#         evald_spec = np.tile(source.eval_spec(lam_arr), ccd.cos_th.size).reshape(ll.T.shape).T
#         lin_out = np.trapz(evald_spec * etalon.eval_airy(ll, cth), lam_arr, axis=0)
#     else:
#         lin_out = np.zeros_like(ccd.cos_th)
#         spec = source.eval_spec(lam_arr)
#         #per_print = (100./ccd.cos_th.size)
#         for i, cth in enumerate(ccd.cos_th):
#             lin_out[i] = np.trapz(spec * etalon.eval_airy(lam_arr, cth), lam_arr)
#             #print "{0:.2f}% done...".format(i * per_print)
#     return lin_out/lin_out.max()
#
# def interpolate_fabry_to_ccd(lin_arr, lin_data, ccd_grid):
#     return griddata(lin_arr, lin_data, ccd_grid, fill_value=0.0)

def eval_fabry(source, ccd, etalon, verbose=False, grid_out=False, both_out=True):
    '''
    Main FP evaluation function. Takes information from input Source class, CCD class, and Etalon class
    to convolve the source function with the Etalon transmission function and then cast it onto the CCD.
    
    :param source: Source class object (example Argon())
    :param ccd: CCD class object (example NikonD5200())
    :param etalon: Etalon class object (example Etalon())
    :param verbose: bool (default=F), if true will print out details during process
    :param grid_out: bool (default=F), if true will return the interpolated image on grid, leave false to avoid interpolation
    :param both_out: bool (default=T), if true and grid_out is true will return both linear array and image
    :return: grid_out=F, lin_data a 1D numpy array of the convolution (think of it as a ring sum) in pixel space
             grid_out=T, both_out=T, lin_data and ccd_data--a 2D numpy array of the ring image in pixel space
             grid_out=T, both_out=F, only the ccd_data
    '''
    lam_arr = np.linspace(source.bounds[0], source.bounds[1], 1000, endpoint=True)
    tic = time.time()
    if ccd.cos_th.size < 1.e5:
        if verbose:
            print 'Evaluating convolution meshgrid...'
        ll, cth = np.meshgrid(lam_arr, ccd.cos_th, indexing='ij')
        if verbose:
            print 'Evaluating source spectrum on meshgrid...'
        evald_spec = np.tile(source.eval_spec(lam_arr), ccd.cos_th.size).reshape(ll.T.shape).T
        if verbose:
            print 'Evaluating convolution integral...'
        lin_out = np.trapz(evald_spec * etalon.eval_airy(ll, cth), lam_arr, axis=0)
    else:
        if verbose:
            print 'CCD linear array has more than 100K points. Using for loop method...'
        lin_out = np.zeros_like(ccd.cos_th)
        spec = source.eval_spec(lam_arr)
        filt = model.eval_spec(lam_arr, 1.0, 468.67, .46 / 2.0 / np.sqrt(2.0 * np.log(2.0)))
        filt *= .3423 / filt.max()
        if True:
            print "Using He filter"
            spec *= filt
        
        #per_print = (100./ccd.cos_th.size)
        for i, cth in enumerate(ccd.cos_th):
            lin_out[i] = np.trapz(spec * etalon.eval_airy(lam_arr, cth), lam_arr)
            #print "{0:.2f}% done...".format(i * per_print)
    lin_data = lin_out / lin_out.max()
    convo_time = time.time() - tic
    if True:
        print "Added gaussian fall off"
        lin_data *= np.exp(-(ccd.lin_arr/3500.0)**2)
    if verbose:
        print 'Convolution complete! It took {0:.2f} seconds.'.format(convo_time)
    if not grid_out:
        return lin_data
    else:
        if verbose:
            print 'Interpolating linear output to ccd grid...'
        tic = time.time()
        ccd_data = griddata(ccd.lin_arr, lin_data, ccd.pix_arr, fill_value=0.0)
        interp_time = time.time() - tic
        if verbose:
            print 'Interpolation complete! It took {0:.2f} seconds.'.format(interp_time)
        if both_out:
            return lin_data, ccd_data
        else:
            return ccd_data

def peak_calc(source, ccd, etalon):
    wav = source.wavelength
    if type(wav) is float:
        m_max = 2.e6 * etalon.d / wav
        num_pks = int(np.floor(np.floor(m_max) - m_max*ccd.cos_th.min() + 1))
        pk_locations = ccd.f * np.sqrt((m_max / (np.floor(m_max) - np.arange(num_pks)))**2 - 1.)
        pk_m = m_max * ccd.f/np.sqrt(ccd.f**2 + pk_locations**2)
    else:
        m_max = np.zeros(len(wav))
        pk_locations = []
        pk_m = []
        for i, w in enumerate(wav):
            mmax = 2.e6 * etalon.d / w
            print w,mmax
            numpks = int(np.floor(np.floor(mmax) - mmax*ccd.cos_th.min() + 1))
            print numpks
            pkloc = ccd.f * np.sqrt((mmax / (np.floor(mmax) - np.arange(numpks)))**2 - 1.)
            print pkloc
            pkm = mmax * ccd.f/np.sqrt(ccd.f**2 + pkloc**2)
            m_max[i] = mmax
            pk_locations.append(pkloc)
            pk_m.append(pkm)
    return {'wavelength': wav, 'm_max': m_max, 'pk_loc': pk_locations, 'pk_m': pk_m}

def main(light='Ar', camera='NikonD5200', amp=None, temp=1.0, vel=0.0, lens=150., d=0.88, F=20., lin_pts=5e4,
         plotit=True, savedic=None, saveccd=None, saveparams=None):

    if light.lower() in ['argon', 'ar']:
        source = Argon(temp=temp, vel=vel)
    elif light.lower() in ['thorium', 'th']:
        source = Thorium()
    elif light.lower() in ['thorium_lamp', 'th_lamp', 'lamp']:
        if amp is None:
            source = Thorium(argon_temp=temp, argon_amp=amp)
        else:
            source = Thorium(argon_temp=temp, argon_amp=amp[0])
    elif light.lower() in ['helium', 'he']:
        source = Helium(vel=vel, temp=temp, rel_amps=amp)
    elif light.lower() in ['mercury', 'hg']:
        source = Mercury(temp=temp, vel=vel)
    else:
        raise Exception('not a valid source name, try again')

    if camera.lower() in ['nikond5200', 'nikon', 'd5200']:
        ccd = NikonD5200(lens=lens, lin_pts=lin_pts)
    elif camera.lower() in ['andor', 'istar', 'istar334t', 'andoristar334t']:
        ccd = Andor(lens=lens, lin_pts=lin_pts)
    else:
        raise Exception('not a valid camera name, try again')

    etalon = Etalon(d=d, F=F)

    pixel_arr = ccd.lin_arr
    costh_arr = ccd.cos_th
    if type(source.wavelength) is float:
        instrument_func = etalon.eval_airy(source.wavelength, costh_arr)
    else:
        instrument_func = etalon.eval_airy(source.wavelength[0], costh_arr)

    pks = peak_calc(source, ccd, etalon)

    lin_data, ccd_data = eval_fabry(source, ccd, etalon, verbose=plotit, grid_out=True)

    param_dict = {"source": source.name, "ccd": ccd.name, "T": temp, "V": vel, "f": (lens, ccd.f), "d": d,
                  "F": F, "Q": etalon.Q, "lam0": source.lam0, "mu": source.mu, "FWHM": source.FWHM,
                  "lam_shfit": source.wavelength, "ccd_type": ccd.name, "npx": ccd.npx, "px_size": ccd.px_size,
                  "ccd_center": ccd.center, "m_max": pks['m_max'], "pk_locations": pks['pk_loc'], "pk_m": pks['pk_m']}

    arr_dict = {"pixel_grid": ccd.pix_arr, "ccd_data": ccd_data, "pixel_arr": ccd.lin_arr, "costh_arr": costh_arr,
                "lin_data": lin_data, "instrument_func": instrument_func}

    if savedic is not None:
        f = h5py.File(savedic, 'w')
        for k in arr_dict.keys():
            f.create_dataset(k, data=arr_dict[k], compression='lzf')
        for k in param_dict.keys():
            f.create_dataset(k, data=param_dict[k])
        print 'Output saved as h5 file: {0}'.format(savedic)

    if saveccd is not None:
        np.save(saveccd, ccd_data)

    if saveparams is not None:
        pickle.dump(param_dict, open(saveparams, 'wb'))

    if plotit:
        print 'Plotting...'
        f, ax = plt.subplots()
        ax.plot(pixel_arr, lin_data, label='data')
        #ax.plot(pixel_arr**2, lin_data, label='data')
        #ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
        #ax.plot(pixel_arr, instrument_func, label='instr. f')
        ax.set_xlabel('pixels')
        #ax.set_xlim([0, min(ccd.npx) / 2.])
        # last_pk_ix = np.where(pk_locations > min(ccd.npx) / 2.)[0][0]
        # ax2 = ax.twiny()
        # ax2.set_xlim(ax.get_xlim())
        # ax2.set_xticks(pk_locations[0:last_pk_ix])
        # ax2.set_xticklabels(pk_m_arr[0:last_pk_ix].astype(str))
        # ax2.set_xlabel(r'm # using $\lambda=${0:.4f}'.format(source.wavelength))
        # for loc in pk_locations:
        #    ax.plot([loc]*2, [0, 1], 'k--')
        ax.legend(loc='best')
        plt.show(block=False)

        f, ax = plt.subplots()
        ax.imshow(ccd_data, cmap='gray', origin='lower')
        ax.set_aspect(1.0)
        plt.show()

        print 'Done!'

    return param_dict, arr_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Forward modeling for Fabry-Perot spectrometer on WiPAL')
    parser.add_argument('light_source', type=str, help='The spectrum that you wish to calculate FP output from.'
                                                    'Current choices are: Argon | Ar, Thorium | Th, '
                                                    'Thorium_lamp | Th_lamp, and Helium | He.')
    parser.add_argument('--camera', '-c', type=str, default='NikonD5200', help='Camera that you wish to use. Current'
                                                                               'choices are Nikon|D5200|NikonD5200 or'
                                                                               'Andor|iStar')
    parser.add_argument('--temp', '-t', type=float, default=1.0, help='Temperature in eV of light soruce. For Thorium'
                                                                      'lamp source, this is the temp of the Argon line')
    parser.add_argument('--vel', '-v', type=float, default=0.0, help='Velocity in km/s of light source. This has no'
                                                                     'effect for the Thorium_lamp source.')
    parser.add_argument('--lens', '-L', type=float, default=150.0, help='Focal length of the camera lens used in mm.')
    parser.add_argument('-d', type=float, default=0.88, help='Etalon spacing in mm.')
    parser.add_argument('--finesse', '-F', type=float, default=20.0, help='Etalon finesse (sharpness of instrument'
                                                                          'function).')
    parser.add_argument('--rel_amp', '-ra', type=float, nargs='+', default=None, help='For the Thorium_lamp case, this'
                                                                                      'is the relative amplitude of'
                                                                                      'the Ar line wrt the Th line. In'
                                                                                      'the Helium case, this is a list'
                                                                                      'of 13 relative amplitudes for'
                                                                                      'the individual lines in the '
                                                                                      '468.6nm complex. The default'
                                                                                      'will be 1 for the Thorium lamp'
                                                                                      'and the amplitudes from the '
                                                                                      'Evans thesis for Helium.')
    parser.add_argument('--pts', '-p', type=float, default=5e4, help='Number of points to use in the array that the '
                                                                       'convolution is calculated with. This number'
                                                                       'should be much larger than the size of the '
                                                                       'sensor array. If larger than 1.e5, then a'
                                                                       'slower for loop method is used to avoid eating'
                                                                       'up too much memory.')
    parser.add_argument('--quiet', '-q', action='store_false', help='Flag to turn off plotting, and command'
                                                                               'line output.')
    parser.add_argument('--save_h5', '-h5', action='store_true', help='Saves a h5 file with all data used'
                                                                                'in this calculation. Use --fname flag'
                                                                                'to have a different filename than the'
                                                                                'default.')
    parser.add_argument('--save_ccd', '-ccd', action='store_true', help='Save a npy file containing just the ccd picture'
                                                                 'output. Use --fname flag to have a different filename'
                                                                 'than the default.')
    parser.add_argument('--save_params', '-params', action='store_true', help='Saves a pickle file of just'
                                                                                         'the parameters from this'
                                                                                         'calc (no arrays).')
    parser.add_argument('--fname', '-f', type=str, default=None, help='Filename for any saves. Default is'
                                                                      'SpecName_pts')
    args = parser.parse_args()
    if args.fname is None:
        if args.save_h5:
            savedic = args.light_source + '_{0:.0f}.h5'.format(args.pts)
        else:
            savedic = None
        if args.save_ccd:
            saveccd = args.light_source + '_{0:.0f}.npy'.format(args.pts)
        else:
            saveccd = None
        if args.save_params:
            saveparams = args.light_source + '_{0:.0f}.pckl'.format(args.pts)
        else:
            saveparams = None
    else:
        if args.save_h5:
            savedic = args.fname + '.h5'
        else:
            savedic = None
        if args.save_ccd:
            saveccd = args.fname + '.npy'
        else:
            saveccd = None
        if args.save_params:
            saveparams = args.fname + '.pckl'
        else:
            saveparams = None

    _, _ = main(light=args.light_source, camera=args.camera, amp=args.rel_amp, temp=args.temp, vel=args.vel,
                lens=args.lens, d=args.d, F=args.finesse, lin_pts=args.pts, plotit=args.quiet, savedic=savedic,
                saveccd=saveccd, saveparams=saveparams)
