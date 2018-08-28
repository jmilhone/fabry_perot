Plasma
=========

Plasma Module
----------------

.. automodule:: fabry.plasma.plasma
    :members:


Argon Plasma Solver Module
------------------------------

.. py:currentmodule:: fabry.plasma.argon_plasma_solver.py

This module contains a few different MultiNest plasma solvers.

.. attribute:: w0

    wavelength for the Ar II ion line

.. attribute:: mu

    mass of argon in amu

.. function:: no_vel_solver(output_folder, Fpost, Lpost, dpost, resume=True, test_plot=False)

    MultiNest solver for argon with no velocity shift

    :param str output_folder: path to folder for input and output folders
    :param np.ndarray Fpost: posterior results for the finesse
    :param np.ndarray Lpost: posterior results for the camera focal length
    :param np.ndarray dpost: posterior results for the etalon spacing
    :param bool resume: resume calculation if True, default=True
    :param bool test_plot: make a test plot of prior intstead of solving, default=False

.. function:: const_vel_solver(output_folder, Fpost, Lpost, dpost, resume=True, test_plot=False)

    MultiNest solver for argon with a constant velocity shift

    :param str output_folder: path to folder for input and output folders
    :param np.ndarray Fpost: posterior results for the finesse
    :param np.ndarray Lpost: posterior results for the camera focal length
    :param np.ndarray dpost: posterior results for the etalon spacing
    :param bool resume: resume calculation if True, default=True
    :param bool test_plot: make a test plot of prior intstead of solving, default=False

.. function:: profile_vel_solver(output_folder, Fpost, Lpost, dpost, resume=True, test_plot=False)

    MultiNest solver for argon with a PCX-U torodial velocity with only outer drive

    :param str output_folder: path to folder for input and output folders
    :param np.ndarray Fpost: posterior results for the finesse
    :param np.ndarray Lpost: posterior results for the camera focal length
    :param np.ndarray dpost: posterior results for the etalon spacing
    :param bool resume: resume calculation if True, default=True
    :param bool test_plot: make a test plot of prior intstead of solving, default=False

