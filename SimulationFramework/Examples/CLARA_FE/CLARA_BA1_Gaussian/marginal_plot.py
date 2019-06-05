import sys, os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
# try:
sys.path.append('C:/Users/jkj62/Documents/GitHub/OnlineModel/')
print sys.path
import SimulationFramework.Modules.read_beam_file as rbf
# except:
#     print 'cant find rbf!'
#     import read_beam_file as rbf
import argparse
parser = argparse.ArgumentParser(description='Analyse genesis output file and return relevant parameters')
parser.add_argument('-d', '--directory', default='.', nargs='+')
parser.add_argument('-g', '--glob', default=None)
parser.add_argument('-f', '--filename', default='CLA-FMS-APER-02.hdf5')
args = parser.parse_args()

def beam_analysis(d):
    beam = rbf.beam()
    print args.filename
    beam.read_HDF5_beam_file(d+'/'+args.filename)
    return beam

def main():
    dirs = args.directory if isinstance(args.directory, (list, tuple)) else [args.directory]

    if args.glob is not None:
        if isinstance(args.directory, (list, tuple)):
            dirs = []
            for d in args.directory:
                dirs += glob(d + '/' + args.glob)
        else:
            dirs = glob(args.directory + '/' + args.glob)

    fig = plt.figure(figsize=(8,8))
    fig.subplots_adjust(wspace=0, hspace=0)
    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_xDist = plt.subplot(gs[0, :2],sharex=ax_main)
    ax_yDist = plt.subplot(gs[1:3, 2],sharey=ax_main)

    for d in dirs:

        beam = beam_analysis(d)
        x, y = 1e12*(beam.t - np.mean(beam.t)), 1e-6*beam.cp
        t = 1e12*(beam.t - np.mean(beam.t))
        dt = 0.05*(max(t) - min(t))

        ax_main.scatter(x,y,marker='.')
        ax_main.set(xlabel="t [ps]", ylabel="cp [MeV/c]")

        beam.slices = 100
        beam.bin_time()
        x = beam.slice['t_Bins']
        bin_centers = 10**12 * np.array((0.5*(x[1:]+x[:-1]) - np.mean(beam.t)))
        ax_xDist.plot(bin_centers, beam.slice_peak_current)
        ax_xDist.set(ylabel='I [A]')

        beam.bin_momentum(width=None)
        y = beam.slice['cp_Bins']
        bin_centers = 10**-6 * np.array((0.5*(y[1:]+y[:-1])))
        # ax_yDist.plot(beam.slice_peak_current, bin_centers)
        normhist = map(lambda x: 100*float(x)/len(beam.t), beam._hist)
        ax_yDist.plot(normhist, bin_centers)
        ax_yDist.set(xlabel='Fraction [%]')

    plt.setp(ax_xDist.get_xticklabels(), visible=False)
    plt.setp(ax_yDist.get_yticklabels(), visible=False)

    plt.show()


if __name__ == '__main__':
    main()
