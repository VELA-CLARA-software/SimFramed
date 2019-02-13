import os, errno, sys
import numpy as np
import random
from scipy.constants import c
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from ocelot.S2E_STFC import FEL_simulation_block
from ocelot.adaptors.genesis import generate_input, get_genesis_launcher, run_genesis, rematch_edist, edist2beam
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
from ocelot.cpbd.beam import Twiss
import SimulationFramework.Modules.read_beam_file as rbf
from SimulationFramework.Modules.optimiser import optimiser
opt = optimiser()
from SimulationFramework.Modules.constraints import constraintsClass
import matplotlib.pyplot as plt
import time
import csv
import multiprocessing
from scoop import futures
from deap import base, creator, tools, algorithms
import copy
import SimulationFramework.Examples.CLARA.Elegant.runElegant as runEle
from shutil import copyfile
from genesisBeamFile import FEL_sim

def read_output_file(dir, run=False):

    # Genesis run:
    data={'gen_file':'NEW_SHORT_NOM_TD_v7.in',
          'file_pout': dir,
          'file_beam':'beamfile.txt',
          'i_scan':0,
          'gen_launch':'',
          'stat_run':1,
          'idump':0,
          'i_edist':0,
          'i_beam':0,
          'i_rad':0,
          'i_dpa':0,
          'i_dfl':0,
          'i_HDF5':1,
          'HDF5_file': dir+'/'+'CLA-FMS-APER-01.hdf5',
          'i_match': 0}
    if run:
        data['gen_launch'] = '/opt/OpenMPI-3.1.3/bin/mpiexec --timeout 300 -np 25 /opt/Genesis/bin/genesis2 < tmp.cmd 2>&1 > /dev/null'
    else:
        data['gen_launch'] = ''
    f = FEL_sim(data)
    A_inp = f.read_GEN_input_file()
    g = f.GEN_simul_preproc(A_inp, dirno=dir)

    p = map(np.mean, zip(*g.p_int))
    g.Lsat, g.Lsatindex = find_saturation(p, g.z)

    # this section is copied from genesis_plot.py - to calculate bandwidth
    spectrum_lamdwidth_fwhm = np.zeros_like(g.z)
    spectrum_lamdwidth_std = np.zeros_like(g.z)
    for zz in range(g.nZ):
        spectrum_lamdwidth_fwhm[zz] = None
        spectrum_lamdwidth_std[zz] = None
        if np.sum(g.spec[:,zz])!=0:
            pos, width, arr = fwhm3(g.spec[:, zz])
            if width != None:
                if arr[0] == arr[-1]:
                    dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos-1])
                else:
                    dlambda = abs( (g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]) )
                spectrum_lamdwidth_fwhm[zz] = dlambda * width / g.freq_lamd[pos]
                # spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / g.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)
            spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz], 0, 1)
    ######### end of section copied from genesis_plot.py
    brightness = g.energy / spectrum_lamdwidth_std

    for i in range(len(g.z)):
        if g.z[i] < 5:
            brightness[i] = 0
    maxb = np.argmax(brightness)
    maxe = np.argmax(g.energy)
    return [ g.energy[maxb], spectrum_lamdwidth_std[maxb], g.z[maxb], g.energy[maxe], spectrum_lamdwidth_std[maxe], g.z[maxe] ]

if __name__ == "__main__":
    with open('reanalyse.txt','w', buffering=0 ) as outfile:
        # csv.register_dialect("custom", delimiter=", ", skipinitialspace=True)
        writer = csv.writer(outfile)
        for i in range(4354):
            tmpdir = './outData/'+str(i)
            sys.stdout = open(tmpdir+'/'+'std.out', 'w')
            sys.stderr = open(tmpdir+'/'+'std.err', 'w')
            ans = read_output_file(tmpdir, run=True)
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            ans = [i] + ans
            print 'ans = ', ans
            writer.writerow(ans)
        print ans
