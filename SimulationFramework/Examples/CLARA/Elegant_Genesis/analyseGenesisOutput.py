import sys, os
sys.path.append(os.path.abspath(__file__+'/../../../../../../S2E_FEL_dev_gen/'))
from ocelot.adaptors.genesis import read_out_file
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
import numpy as np
import argparse
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf

parser = argparse.ArgumentParser(description='Analyse genesis output file and return relevant parameters')
parser.add_argument('-d', '--directory', default='.')
parser.add_argument('-f', '--file', default='run.0.gout')

def analysis(g):
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
    brightness[0] = 0
    return g.energy, spectrum_lamdwidth_std, g.z, brightness

def beam_analysis(d):
    beam = rbf.beam()
    beam.read_HDF5_beam_file(d+'/CLA-FMS-APER-01.hdf5')
    return beam

def main():
    args = parser.parse_args()
    out = read_out_file(args.directory + '/' + args.file, read_level=2)
    e, b, l, bright = analysis(out)
    mb = np.argmax(bright)
    print 'bandwidth[br] = ', 1e2*b[mb], '%  pulse energy[br] =', 1e6*e[mb], 'uJ  Sat. Length[br] =', l[mb], 'm  Brightness[max] = ', bright[mb]
    en = np.argmax(e)
    print 'bandwidth[en] = ', 1e2*b[en], '%  pulse energy[max] =', 1e6*e[en], 'uJ  Sat. Length[en] =', l[en], 'm  Brightness[en] = ', bright[en]
    beam = beam_analysis(args.directory)
    print 'bunchlength = ', 1e6*beam.sigma_z, 'um  chirp = ', -1*beam.chirp, 'MeV/ps'


if __name__ == '__main__':
   main()
