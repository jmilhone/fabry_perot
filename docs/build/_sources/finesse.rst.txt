Finesse
==========

Argon Solver
-------------

.. py:currentmodule:: fabry.finesse.argon_solver

This module contains a few different MultiNest calibration solvers for the Th lamp with the 488 nm central wavelength 1 nm fwhm filter. 

.. attribute:: mu  

    Tuple with the relative masses for Th and  Ar

.. attribute:: w0

    Tuple with Th I and Ar II wavelengths used for calibration

.. function:: solver(output_folder,prior_filename,data_filename,Lpost,dpost,resume=True,test_plot=False)

    MultiNest solver for point spread function calibration with the Argon filter

    :param str output_folder: path to folder containing input and output files for finesse solver
    :param str prior_filename: path to file for the MultiNest prior (json format)
    :param str data_filename: path to file for the input data 
    :param np.ndarray Lpost: array of equally weighted marginal posterior L values
    :param np.ndarray dpost: array of equally weighted marginal posterior d values
    :param resume (bool, optional): MultiNest will resume the calculation if True
    :param test_plot: for debugging purposes, allows the user to sample the prior and compare to input data
    :type test_plot: bool or None
    :return: None

.. function:: full_solver(output_folder,prior_filename,data_filename,resume=True,test_plot=False)

    MultiNest solver for point spread function calibration with the Argon filter. This is a 
        full solver. L and d will be solved as well!

    :param str output_folder: path to folder containing input and output files for finesse solver
    :param str prior_filename: path to file for the MultiNest prior (json format)
    :param str data_filename: path to file for the input data 
    :param resume: MultiNest will resume the calculation if True
    :type resume: bool or None
    :param test_plot: for debugging purposes, allows the user to sample the prior and compare to input data
    :type test_plot: bool or None


Check Argon Solver
--------------------

.. py:currentmodule:: fabry.finesse.check_argon_solver

.. attribute:: mu  

    Tuple with the relative masses for Th and  Ar

.. attribute:: w0

    Tuple with Th I and Ar II wavelengths used for calibration


.. function:: check_solver(finesse_folder, Lpost, dpost)

    Verfies the output of the point spread function solver.  Makes many plots.

    :param str finesse_folder: path to folder containing output files from the finesse solver
    :param np.ndarray Lpost: array of equally weighted marginal posterior L values
    :param np.ndarray dpost: array of equally weighted marginal posterior d values

.. function:: check_full_solver(finesse_folder)

    Verfies the output of the full point spread function solver.  Makes many plots.

    :param str finesse_folder: path to folder containing output files from the finesse solver



