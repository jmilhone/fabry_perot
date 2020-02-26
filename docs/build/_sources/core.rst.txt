Core
=======

RingSum
----------

.. automodule:: fabry.core.ringsum
    :members: get_bin_edges, locate_center, ringsum

.. py:currentmodule:: fabry.core.ringsum

.. function:: super_pixelate(data, npix=2)

    Creates super pixels for image data

    :param np.ndarray data: 2d image data
    :param int npix: integer number of pixels to create an npix x npix super pixel, default=2
    :returns: new image made from the super pixels
    :rtype: np.ndarray

.. function:: calculate_weighted_mean(data, error)

    Calculates the weighted mean of data with standard deviation error

    :param np.ndarray data: some stuff
    :param np.ndarray error: standard deviation error bar for data
    :return: weighted mean and weighted standard deviation
    :rtype: tuple (float, float)

Models
----------

.. automodule:: fabry.core.models
    :members:

Fitting
------------

.. automodule:: fabry.core.fitting
    :members:

