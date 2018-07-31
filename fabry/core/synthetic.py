from __future__ import print_function, division, absolute_import
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import argparse
from . import models
import multiprocessing as multi
from collections import MutableMapping
import numbers


########################################
# Physical constants, DO NOT OVERWRITE #
########################################
q = 1.6e-19
c = 2.998e8
mp = 1.672e-27


class Sensor(object):
    """A representation for a camera sensor for the Fabry Perot

    Attributes:
        nx (numbers.Integral): number of pixels in the x direction
        ny (numbers.Integral): number of pixels in the y direction
        px_size (numbers.Real): Pixel size in mm
        x0 (numbers.Real): x location of the center
        y0 (numbers.Real): y location of th center
        sensor (np.ndarray): 2D array representing the sensors where each 
            element is the count value for that pixel
    """
    def __init__(self, nx, ny, px_size=0.004):
        super(Sensor, self).__init__()

        self.nx = int(nx)
        self.ny = int(ny)
        self.sensor = np.zeros((nx, ny))
        self.px_size = px_size
        self.x0 = nx/2.0
        self.y0 = ny/2.0
        self._R = None  # Delay the creation of the R matrix until it is actually needed

    def create_Rgrid(self):
        """Creates a 2D matrix of radii values for each pixel

        returns:
            np.ndarray
        """
        nx = self.nx
        ny = self.ny

        x = np.arange(1, nx+1, 1)
        y = np.arange(1, ny+1, 1)

        XX, YY = np.meshgrid(x, y)

        R = np.sqrt( (XX-self.x0)**2 + (YY-self.y0)**2 )

        return R

    @property
    def R(self):
        """np.ndarray: 2D array of radii values for each pixel"""
        if self._R is None:
            self._R = self.create_Rgrid()
        return self._R

    @R.setter
    def R(self):
        self._R = self.create_Rgrid()

    def calculate_emission(self, etalon, light_source, nprocs=4):
        """Calculates emission from a light source through an etalon onto the sensor

        Args:
            etalon (Etalon): representation of the etalon
            light_source (LightSource): representation of the light source
        """
        self.sensor = etalon.calculate_emission(self, light_source, nprocs=nprocs)

    @classmethod
    def from_dict(cls, sensor):
        """Creates an instance of Sensor from a dictionary

        Args:
            sensor (dict): dictionary representing a sensor

        Returns:
            Sensor: a new instance of a Sensor
        """
        nx = sensor.get('nx', 1024)
        ny = sensor.get('ny', 1024)
        px_size = sensor.get('px_size', 0.004)
        return cls(nx, ny, px_size=px_size)

    def to_dict(self):
        """Returns a dictionary representation of a Sensor

        Returns:
            dict: a dictionary representation of a Sensor
        """
        return {'nx': self.nx, 'ny': self.ny, 'px_size': self.px_size}

    def __repr__(self):
        class_name = type(self).__name__
        return '{}({!r}, {!r}, px_size={!r})'.format(class_name, self.nx, self.ny, self.px_size)


class Etalon(object):
    """
    Class that represents an etalon for a Fabry-Perot spectrometer

    Attributes:
        L (numbers.Real): focal length of lens for the camera
        d (numbers.Real): etalon spacing
        F (numbers.Real): finesse of etalon

    """
    def __init__(self, L, d, F):
        super(Etalon, self).__init__()
        self.L = L
        self.d = d
        self.F = F

    def calculate_emission(self, sensor, light_source, nprocs=4):
        """Calcultes the emission onto a sensor from a light source

        Note: This uses multiprocessing for speed. Each process has a for loop 
            because of memory constraints if the sensor is too big.

        Args:
            sensor (Sensor): represents a camera sensor
            light_source (LightSource): represents a light source for the Fabry Perot

        kwargs:
            nprocs (numbers.Integral): number of processes to use

        Returns:
            np.ndarray: shape matches sensor.sensor.shape
        """
        r = sensor.R
        px_size = sensor.px_size

        w = light_source.wavelength
        mu = light_source.mu
        amp = light_source.amplitude
        vel = light_source.velocity
        temp = light_source.temperature

        split_r = np.array_split(r.flatten(), nprocs)
        procs = []
        sigs = {}
        out = multi.Queue()
        labels = ['{0}'.format(x) for x in range(nprocs)]
        for k in range(nprocs):
            p = multi.Process(target = Etalon._calculate_emission,
                    args=(split_r[k], self.L / px_size, self.d, self.F, w, mu, amp, temp, vel),
                    kwargs={'out': out, 'label': labels[k]})
            procs.append(p)
            p.start()

        for i in range(nprocs):
            tup = out.get()
            sigs[tup[0]] = tup[1]

        for p in procs:
            p.join()

        emission = []
        for k in labels:
            emission.append(sigs[k])
        emission = np.concatenate(emission)
        emission.shape = r.shape

        return emission

    @staticmethod
    def _calculate_emission(r, L, d, F, w, mu, amp, temp, vel, out=None, label=None, noise=True):
        """Helper function for calculating emission

        Note: This is utlized by calculate_emission for multiprocessing

        Args:
            r (np.ndarray): radii in pixels
            L (numbers.Real): focal length of camera lens in pixels
            d (numbers.Real): etalon spacing in mm
            F (numbers.Real): finesse of etalon
            w (numbers.Real): wavelength in nm
            mu (numbers.Real): relative mass
            amp (numbers.Real): amplitude of line
            temp (numbers.Real): ion temperature in eV
            vel (numbers.Real): velocity of ion in m/s

        Kwargs:
            out (multiprocessing.Queue): output queue
            label (str): label for the output being put into the output queue

        Returns:
            np.ndarray (optional) only if not using for multiprocessing
        """
        # Memory is a problem here because of broadcasting so I'm going to split the problem up
        model = []
        r_split = np.array_split(r.flatten(), 1000)
        for r_sp in r_split:
            model.append(models.forward_model(r_sp, L, d, F, w, mu, amp, temp, vel))
        model = np.concatenate(model)

        if noise:
            print(label, 'adding noise to the image')
            npts = len(model)
            for i in xrange(npts):
                model[i] = model[i] + np.random.normal(scale=np.sqrt(model[i]))
                if i % 100 ==0:
                    print(model[i], np.random.normal(scale=np.sqrt(model[i])))

        if out and label:
            print(label)
            out.put((label, model))
        else:
            return model

    def __repr__(self):
        class_name = type(self).__name__
        return "{}({!r}, {!r}, {!r})".format(class_name, self.L, self.d, self.F)

    @classmethod
    def from_dict(cls, etalon):
        """Creates an instance of Etalon from a dictionary

        Args:
            etalon (dict): dictionary representing a etalon

        Returns:
            Etalon: an new instance of a Etalon
        """
        L = etalon.get('L', 150.0)
        d = etalon.get('d', 0.88)
        F = etalon.get('F', 21.0)
        return cls(L, d, F)

    def to_dict(self):
        """Returns a dictionary representing itself

        Returns:
            dict: a dictionary representing itself
        """
        return {"L": self.L, "d": self.d, "F": self.F}


class LightSource(object):
    """A representation of a light source for a Fabry-Perot spectrometer

    Attributes:
        temperature (numbers.Real): temperature of the emitting ion in eV
        wavelength (numbers.Real): wavelength of the light emitted in nm
        mu (numbers.Real): relative mass of the ion
        amplitude (numbers.Real): amplitude of the light emitted (you can choose your units here...)
        velocity (VelocityProfile or numbers.Real): velocity of the emitting ion in m/s
    """

    def __init__(self, Ti, w, mu, velocity, amplitude=1):
        super(LightSource, self).__init__()
        self.temperature =  Ti
        self.wavelength = w
        self.mu = mu
        self.amplitude = amplitude
        self.velocity = velocity

    def __repr__(self):
        class_name = type(self).__name__
        return "{}({!r}, {!r}, {!r}, {!r}, amplitude={!r})".format(
                class_name, self.temperature, self.wavelength, self.mu, 
                self.velocity, self.amplitude)

    @classmethod
    def from_dict(cls, light_source):
        """Creates a new instance of LightSource from a dictionary

        Args:
            light_source (dict): dictionary representing a light source

        Returns:
            LightSource
        """
        temperature = light_source.get('temperature', 0.5)
        wavelength = light_source.get('wavelength', 488.0)
        mu = light_source.get('mu', 40.0)
        amplitude = light_source.get('amplitude', 1.0)

        # Oomph, not sure I like this, but I couldn't think of a better way
        velocity = light_source.get('velocity', 0.0)
        if isinstance(velocity, MutableMapping):
            vel_class = velocity.get("class_name", VelocityProfile)
            vel_class = globals().get(vel_class, None)
            velocity.pop('class_name')
            if vel_class is None:
                velocity = 0.0
            else:
                velocity = vel_class(**velocity)

        return cls(temperature, wavelength, mu, velocity, amplitude=amplitude)

    def to_dict(self):
        """Returns a dict representing itself"""
        velocity = 0.0
        if self.velocity is not None:
            velocity = self.velocity

        try:
            velocity = velocity.to_dict()
        except AttributeError:
            pass

        return {
                'temperature': self.temperature,
                'wavelength': self.wavelength,
                'mu': self.mu,
                'amplitude': self.amplitude,
                'velocity': velocity,
                }



class UniformPlasma(LightSource):
    """A representation of a uniform density and uniform Ti plasma with ion 
        species mu emitting light for a Fabry-Perot spectrometer

    Attributes:
        temperature (numbers.Real): temperature of the emitting ion in eV
        wavelength (numbers.Real): wavelength of the light emitted in nm
        mu (numbers.Real): relative mass of the ion
        amplitude (numbers.Real): amplitude of the light emitted (you can choose your units here...)
        velocity (Union[VelocityProfile,numbers.Real]): velocity of the emitting ion in m/s
        ne (numbers.Real): electron density in cm^-3
        pec (numbers.Real): photon emissivity coefficient (need to decide on units here)
        mu (numbers.Real): relative mass of the ion
    """

    def __init__(self, ne, Ti, pec, w, velocity=None, mu=40.0):
        super(UniformPlasma, self).__init__(Ti, w, mu, velocity)

        self.ne = ne
        self.pec = pec
        self.mu = mu

    def ion_emission(self, r, wavelength, cos_theta=None):
        """Calculates ion emission at a radius r for wavelengths provided

        Args:
            r (Union[numbers.Real, np.ndarray]): radii to calculate ion emission at
            wavelength (Union[numbers.Real, np.ndarray]): wavelength to calculate emission line 
                profile

        Kwargs:
            cos_theta (Union[numbers.Real, np.ndarray]): cos(theta) to project velocity onto a unit
                vector an angle theta from the toroidal direction

        Returns:
            np.ndarray
        """
        velocity = 0.0
        if self.velocity is not None:
            velocity = self.velocity

        if callable(velocity):
            velocity = velocity(r)

        if cos_theta is not None:
            print('im in the cos theta portion')
            cosine = np.asarray(cos_theta)
            print(cosine.max(), cosine.min())
            velocity *= cosine

        line_profile = self.gaussian(wavelength, velocity)

        emission = self.ne**2 * self.pec * line_profile / (4*np.pi)
        return emission

    def chord_emission(self, impact_factor, wavelength):

        b = np.asarray(impact_factor)
        if np.max(b) < 0.0:
            raise ValueError('impact_factor must be greater than or equal to zero')

        max_radii = 150.0
        x_max = np.sqrt(max_radii**2 - b**2)
        x_arr = np.linspace(0.0, x_max, 1000)

        # I need the x_arr and subsequent arrays to be broadcastable with wavelength array
        x_arr = x_arr[np.newaxis, :]
        w = wavelength[:, np.newaxis]

        # theta_arr = np.arctan2(b, x_arr)
        cos_theta = b / np.sqrt(b**2 + x_arr**2)

        rarr = np.sqrt(x_arr**2 + b**2)
        print(rarr.shape)
        print(w.shape)
        emission = self.ion_emission(rarr, w, cos_theta=cos_theta)

        radiance = 2.0*np.trapz(emission, x=x_arr, axis=1)
        print(emission.shape, x_arr.shape, w.shape, radiance.shape)
        #fig, ax = plt.subplots()
        #for i in range(1000):
        #    if i % 50 == 0:
        #        ax.plot(wavelength, emission[:, i] / emission.max())
        #ax.plot(wavelength, radiance / radiance.max(), 'k')
        #plt.show()
        return radiance


    def gaussian(self, wavelength, velocity):
        """Calculates doppler broadened and shifted gaussian

        Args:
            wavelength (Union[numbers.Real, np.ndarray]): wavelength to calculate emission line profile
            velocity (Union[numbers.Real, np.ndarray]): velocity of ion for doppler shift
        """
        w = np.asarray(wavelength)
        v = np.asarray(velocity)
        sigma = self.sigma
        w_shift = self.calculate_shift(v)
        norm = np.sqrt(2*np.pi) * sigma

        return np.exp(-0.5*(w-w_shift)**2 / sigma**2) / norm

    @property
    def sigma(self):
        """Thermal doppler broadening"""
        return np.sqrt(q * self.temperature / (self.mass)) * self.wavelength / c

    def calculate_shift(self, velocity):
        """Calculate doppler shift from the ion velocity

        Args:
            velocity (Union[numbers.Real, np.ndarray]): velocity in m/s

        Returns:
            np.ndarray
        """
        return self.wavelength * (1.0 - velocity / c)

    @property
    def mass(self):
        """Mass of ion in kg"""
        return self.mu * mp

    def __repr__(self):
        class_name = type(self).__name__
        return "{}({!r}, {!r}, {!r}, {!r}, velocity={!r}, mu={!r})".format(
                class_name, self.ne, self.temperature, self.pec, self.wavelength,
                self.velocity, self.mu)


    def to_dict(self):
        """Returns a dict representation
        """
        velocity = 0.0
        if self.velocity is not None:
            velocity = self.velocity

        try:
            velocity = velocity.to_dict()
        except AttributeError:
            pass

        return {
                'temperature': self.temperature,
                'wavelength': self.wavelength,
                'mu': self.mu,
                'velocity': velocity,
                'pec': self.pec,
                'ne': self.ne
                }



    @classmethod
    def from_dict(cls, plasma):
        """Creates a new instance of UniformPlasma from dict

        Args:
            plasma (dict): dictionary representation of a UniformPlasma

        Returns:
            UniformPlasma
        """
        temperature = plasma.get('temperature', 0.5)
        wavelength = plasma.get('wavelength', 488.0)
        mu = plasma.get('mu', 40.0)
        pec = plasma.get('pec')
        ne = plasma.get("ne", 1e12)

        # Oomph, not sure I like this, but I couldn't think of a better way
        velocity = plasma.get('velocity', 0.0)
        if isinstance(velocity, MutableMapping):
            vel_class = velocity.get("class_name", VelocityProfile)
            vel_class = globals().get(vel_class, None)
            velocity.pop('class_name')
            if vel_class is None:
                velocity = 0.0
            else:
                velocity = vel_class(**velocity)

        return cls(ne, temperature, pec, wavelength, velocity=velocity, mu=mu)



class VelocityProfile(object):
    """Represents a edge driven velocity profile

    Attributes:
        Vmax (numbers.Real): maximum velocity
        max_radius (numbers.Real): radial location of the maximum velocity
        length_scale (numbers.Real): scale for the velocity gradient inward radially
        edge_scale (numbers.Real): scale for the velocity edge gradient
        R0 (numbers.Real): location of the edge
        offset (numbers.Real): offset velocity in the center
    """
    def __init__(self, Vmax, max_radius, length_scale, edge_scale, R0=140.0, offset=0.0):
        super(VelocityProfile, self).__init__()

        self.Vmax = Vmax
        self.max_radius = max_radius
        self.length_scale = length_scale
        self.edge_scale = edge_scale
        self.R0 = R0
        self.offset=offset

    def vphi(self, r):
        """Returns the Torodial velocity at r

        Args:
            r (Union[numbers.Real, np.ndarray]): radii to evaluate vphi at

        Returns:
            np.ndarray
        """
        radii = np.asarray(r)
        right_profile = self.right_edge_profile(radii)
        left_profile = self.gaussian(radii)

        vel = right_profile * left_profile
        vel *= self.Vmax / np.max(vel)

        return vel

    def right_edge_profile(self, r):
        """Helper function for the edge profile

        Args:
            r (Union[numbers.Real, np.ndarray]): radii to evaulate at

        Returns:
            np.ndarray
        """
        return 0.5 * (1.0 - np.tanh((r - self.R0)/self.edge_scale))

    def gaussian(self, r):
        """Helper function for the inward velocity gradient

        Args:
            r (Union[numbers.Real], np.ndarray]): radii to evaluate at

        Returns:
            np.ndarray
        """
        g = (1-self.offset / self.Vmax)*np.exp(-(r - self.max_radius)**2 / self.length_scale**2)
        g += self.offset / self.Vmax
        #print(self.offset / self.Vmax, 1 + self.offset / self.Vmax)
        #fig, ax = plt.subplots()
        #ax.plot(r, g)
        #plt.show()
        return g

    def __call__(self, r):
        return self.vphi(r)

    def __repr__(self):
        cls = type(self).__name__
        s = "{}({!r},{!r},{!r},{!r},R0={!r},offset={!r})".format(cls, self.Vmax, self.max_radius, self.length_scale, self.edge_scale, self.R0, self.offset)
        return s

    def to_dict(self):
        """Returns a dict representation of the VelocityProfile

        Returns:
            dict
        """

        output_dict = {'class_name': type(self).__name__,
                       'Vmax': self.Vmax,
                       'max_radius': self.max_radius,
                       'length_scale': self.length_scale,
                       'edge_scale': self.edge_scale,
                       'R0': self.R0,
                       'offset': self.offset,
                       }
        return output_dict

    @classmethod
    def from_dict(cls, velocity):
        """Returns a new instance of VelocityProfile from a dict

        Args:
            velocity (dict): dict reprsentation of a VelocityProfile

        Returns:
            VelocityProfile
        """
        Vmax = velocity.get('Vmax')
        max_radius = velocity.get('max_radius')
        length_scale = velocity.get('length_scale')
        edge_scale = velocity.get('edge_scale')
        R0 = velocity.get('R0', 140.0)
        offset = velocity.get('offset', 0.0)
        return cls(Vmax, max_radius, length_scale, edge_scale, R0=R0, offset=offset)


