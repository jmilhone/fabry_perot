from __future__ import print_function, division, absolute_import
import pymultinest
import numpy
from ..core.models import forward_model
import json
import argparse
import h5py
import matplotlib.pyplot as plt # For testing purposes only
from ..tools import file_io as io


def solver(prior_filename, data_filename):

    def log_prior(cube, ndim, nparams):
        pass

    def log_likelihood(cube, ndim, nparams):
        pass

    with open(prior_filename, 'r') as infile:
        prior = json.load(infile, parse_float=np.float64)

    data = io.h5_2_dict(data_filename)


if __name__ == "__main__":
    pass












