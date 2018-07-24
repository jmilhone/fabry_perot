from __future__ import print_function, division
import os
import numpy as np
import h5py

def dict_2_h5(fname, dic, append=False):
    '''
    Writes a dictionary to a hdf5 file with given filename
    It will use lzf compression for all numpy arrays
    
    Args:
        fname (str): filename to write to
        dic (dictionary): dictionary to write
        append (bool, default=False): if true, will
            append to file instead of overwriting
    '''
    if append:
        method = 'r+'
    else:
        method = 'w'
    with h5py.File(fname, method) as h5:
        recursive_save_dict_to_h5(h5, '/', dic)

def h5_2_dict(fname):
    '''
    Reads a dictionary from a hdf5 file with given filename

    Args:
        fname (str): hdf5 filename to read
    Returns:
        dic (dictionary): dictionary of hdf5 keys
    '''
    with h5py.File(fname, 'r') as h5:
        return recursive_load_dict_from_h5(h5, '/')

def prep_folder(path):
    '''
    Checks if folder exists and recursively creates folders
    to ensure the path is valid

    Args:
        path (str): path to folder
    '''
    if os.path.isdir(path):
        return
    else:
        os.makedirs(path)

def recursive_save_dict_to_h5(h5, path, dic):
    ''' function used in save_dict_to_h5 in order to get recursion
    '''
    for key, item in dic.items():
        if path+key in h5:  ### overwrites pre-existing keys with same name
           del h5[path+key]
        if type(item) in [np.ndarray, np.generic]:
            h5.create_dataset(path+key, data=item, compression='lzf')
        elif type(item) != dict:
            try:
                h5.create_dataset(path+key, data=item)
            except TypeError:
                raise ValueError('Cannot save %s type'%type(item))
        else:
            recursive_save_dict_to_h5(h5, path+key+'/', item)

def recursive_load_dict_from_h5(h5, path):
    ''' function used in load_h5_to_dict in order to get recursion
    '''
    out_dict = {}
    for key, item in h5[path].items():
        if type(item) == h5py._hl.dataset.Dataset:
            out_dict[key] = item.value
        elif type(item) == h5py._hl.group.Group:
            out_dict[key] = recursive_load_dict_from_h5(h5, path+key+'/')
    return out_dict

def read_Ld_results(Ld_directory):
    '''
    reads L and d histogram data from multinest run
    
    Args:
        Ld_directory (str): path to multinest save directory
    Returns:
        L (np.ndarray): L histogram values (in pixels)
        d (np.ndarray): d histogram values (in mm)
    '''
    try:
        fname = os.path.join(Ld_directory,"Ld_post_equal_weights.dat")
        post = np.loadtxt(fname, ndmin=2)
    except IOError:
        fname = os.path.join(Ld_directory,"Ld_solver_post_equal_weights.dat")
        post = np.loadtxt(fname, ndmin=2)

    L = post[:, 0]
    d = post[:, 1]
    return L, d

def read_match_finesse_results(finesse_directory, errtemp=False):
    fname = os.path.join(finesse_directory,"F_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    F = post[:, 0]
    V = post[:, 1]
    T = post[:, 2]
    if errtemp:
        E = post[:, 3]
        return F,V,T,E
    else:
        return F,V,T


def read_finesse_results(finesse_directory):
    fname = os.path.join(finesse_directory, "finesse_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    F = post[:, 0]
    A = post[:, 1]
    Arel = post[:, 2]
    Ti = post[:, 3]

    return F, A, Arel, Ti


def read_lyon_temp_results(temp_directory):
    fname = os.path.join(temp_directory,'temp_post_equal_weights.dat')
    post = np.loadtxt(fname, ndmin=2)
    T = post[:,0]
    V = post[:,1]
    #A = post[:,2]
    #O = post[:,3]
    return T,V#,A#,O
