import numpy as np
import matplotlib.pyplot as plt
from os.path import join
import argparse
import sys
sys.path.append('../')
from tools.images import get_data
from tools.file_io import h5_2_dict, dict_2_h5, prep_folder
from core.fitting import find_maximum
from core.models import lyon_temp_forward
import pymultinest
from mpi4py import MPI
from match_finesse_solver import temp_match_solver, match_finesse_check

def write_slice(folder, center=None, atol=1.e-3, sm_range=[515.,770.]):
    a = h5_2_dict(join(folder,'ringsum.h5'))
    data = get_data(a['fname'],color=a['color'])

    R,Z = np.meshgrid(np.arange(data.shape[1],dtype='float64'), np.arange(data.shape[0], dtype='float64'))
    if center is None:
        center = a['center']
    R -= center[0]
    Z -= center[1]

    Theta = -1.*(np.arctan2(Z,-R)-np.pi)
    Rho = np.sqrt(R**2+Z**2)
    del R,Z

    Rho = Rho.flatten()
    Theta = Theta.flatten()
    data = data.flatten()
    ixs = np.where(np.logical_and( Rho < sm_range[1], Rho > sm_range[0]))
    Rho = Rho[ixs]
    Theta = Theta[ixs]
    data = data[ixs]
    thetas = np.arange(0.,360.)
    tsz = thetas.size

    fig,axs = plt.subplots(figsize=(16,9))
    out_data = {'thetas': thetas, 'center':center}
    for i,t in enumerate(thetas):
        print '{0} of {1}'.format(i+1, tsz)
        ixs = np.isclose(Theta, np.ones_like(Theta)*np.deg2rad(t), atol=atol)
        rr = Rho[ixs]
        dd = data[ixs]
        ix = np.argsort(rr)
        axs.plot(rr[ix],dd[ix])
        out_data['{0:0d}'.format(int(t))] = dd[ix]
        out_data['{0:0d}_r'.format(int(t))] = rr[ix]

    dict_2_h5(join(folder,'ring_slices.h5'), out_data)
    plt.show()

def run_Ti_solver():
    #write_slice('../Data/2017_12_11/0199/')
    folder = '../Data/2017_12_11/0199/'
    savefolder = join(folder,'slice_temps/')
    prep_folder(savefolder)
    data = h5_2_dict(join(folder,'ring_slices.h5'))
    Ld_folder = '../Data/2017_12_11/calib/0176'
    F_folder = '../Data/2017_12_11/calib/0179'
    Lpost, dpost = read_Ld_results(Ld_folder)
    Fpost = read_F_results(F_folder)
    current = 800
    
    for i,t in enumerate(data['thetas']):
        if i > 0:
            sig = np.array(data['{0:0d}'.format(int(t))])
            r = np.array(data['{0:0d}_r'.format(int(t))])
            error = 0.05*np.mean(sig) * np.ones_like(sig) + 0.05*sig
            basename = join(savefolder,'{0:0d}_'.format(int(t)))
            temp_solver(r, sig, error, current, Lpost, dpost, Fpost, basename)

def read_Ld_results(Ld_directory):
    '''
    reads L and d histogram data from multinest run
    
    Args:
        Ld_directory (str): path to multinest save directory
    Returns:
        L (np.ndarray): L histogram values (in pixels)
        d (np.ndarray): d histogram values (in mm)
    '''
    fname = join(Ld_directory,"Ld_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    L = post[:, 0]
    d = post[:, 1]
    return L, d

def check_Ti_results():
    folder = '../Data/2017_12_11/0199/'
    data = h5_2_dict(join(folder,'ring_slices.h5'))
    thetas = data['thetas'][1::]
    del data
    savefolder = join(folder,'slice_temps/')
    temp = np.zeros_like(thetas)
    temp_sd = np.zeros_like(thetas)
    for i,t in enumerate(thetas):
        fname = join(savefolder,'{0:0d}_post_equal_weights.dat'.format(int(t)))
        post = np.loadtxt(fname, ndmin=2)
        temp[i] = np.mean(post[:, 0])
        temp_sd[i] = np.std(post[:,0])
    f,axs = plt.subplots(figsize=(16,9))
    axs.errorbar(thetas, temp, yerr=temp_sd, fmt='o')
    plt.show()

def run_match_solver():
#    write_slice('../Data/2017_12_11/0179/')
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        folder = '../Data/2017_12_11/0179/'
        savefolder = join(folder, 'slice_finesse/')
        prep_folder(savefolder)
        data = h5_2_dict(join(folder, 'ring_slices.h5'))
        Ld_folder = '../Data/2017_12_11/calib/0176'
        Lpost, dpost = read_Ld_results(Ld_folder)

        #f,axs = plt.subplots()
        #axs.plot(data['180_r'],data['180'],'o')
        #plt.show()

        sig = np.array(data['180'])
        r = np.array(data['180_r'])
        error = 0.05*np.mean(sig) * np.ones_like(sig) + 0.05*sig
        basename = join(savefolder, 'match_finesse_')
        inputs = {'r':r, 'sig':sig, 'error':error, 'error_per':0.05, 'F_lim':[1.,10.],
                'temp':0.15, 'temp_sigma':0.015, 'wavelength':487.98634, 'mu':39.948, 'livepoints':500,
                'V_lim':[-2.,2.], 'Ld_folder':Ld_folder, 'errtemp':False}
        dict_2_h5(join(savefolder,'input_match_finesse.h5'),inputs)

        solver_in = {'r':r, 'sig':sig, 'error':error, 'Lpost':Lpost, 'dpost':dpost, 'basename':basename}
    else:
        solver_in = None
    
    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        temp_match_solver(solver_in['r'], solver_in['sig'], solver_in['error'],0.15,0.015,solver_in['Lpost'], solver_in['dpost'], solver_in['basename'], F_lim=[1.,10.], errtemp=False) 

    if rank == 0:
        match_finesse_check(savefolder)

if __name__ == "__main__":
    run_match_solver()

