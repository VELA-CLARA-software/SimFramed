import sys, os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
try:
    sys.path.append('C:/Users/jkj62/Documents/GitHub/OnlineModel/')
    import SimulationFramework.Modules.read_beam_file as rbf
    import SimulationFramework.Modules.read_twiss_file as rtf
except:
    import read_beam_file as rbf
import argparse
parser = argparse.ArgumentParser(description='Analyse genesis output file and return relevant parameters')
parser.add_argument('-f', '--file', default='./CLA-FMS-APER-02.hdf5')
parser.add_argument('-d', '--directory', default='.', nargs='+')

def beam_analysis(d, filename):
    beam = rbf.beam()
    beam.read_HDF5_beam_file(d+'/'+filename)
    return beam

def main():
    args = parser.parse_args()
    dirs = args.directory if isinstance(args.directory, (list, tuple)) else [args.directory]

    twiss = rtf.twiss()

    # fig = plt.figure(figsize=(15,15))
    # fig.subplots_adjust(wspace=0, hspace=0)
    # gs = gridspec.GridSpec(3, 3)
    # ax_main = plt.subplot(gs[1:3, :2])
    # ax_xDist = plt.subplot(gs[0, :2],sharex=ax_main)
    # ax_yDist = plt.subplot(gs[1:3, 2],sharey=ax_main)

    for d in dirs:

        beam = beam_analysis(d, args.file)
        twiss.read_elegant_twiss_files(d+'/CLARAX.twi' )
        x, y = 1e12*(beam.t - np.mean(beam.t)), 1e-6*beam.cp
        t = 1e12*(beam.t - np.mean(beam.t))
        dt = 0.05*(max(t) - min(t))

        # ax_main.scatter(x,y,marker='.')
        # ax_main.set(xlabel="t [ps]", ylabel="cp [MeV/c]")

        beam.slices = 100
        beam.bin_time()
        x = beam.slice['t_Bins']
        bin_centers = 10**12 * np.array((0.5*(x[1:]+x[:-1]) - np.mean(beam.t)))
        # ax_xDist.plot(bin_centers, beam.slice_peak_current)
        # ax_xDist.set(ylabel='I [A]')

        fig2, ax = plt.subplots(3, 2,
                               figsize=(10, 10))
        fig2.subplots_adjust(hspace=0.01)

        t = 1e12*(beam.t-np.mean(beam.t))
        t_grid = np.linspace(min(t), max(t), 2**8)
        peakIPDF = beam.PDF(t, t_grid, bandwidth=beam.rms(t)/(2**4))
        peakICDF = beam.CDF(t, t_grid, bandwidth=beam.rms(t)/(2**4))
        peakIFWHM, indexes = beam.FWHM(t_grid, peakIPDF, frac=0.1)
        peakIFWHM2, indexes2 = beam.FWHM(t_grid, peakIPDF, frac=20)
        stdpeakIPDF = (max(peakIPDF[indexes2]) - min(peakIPDF[indexes2]))/np.mean(peakIPDF[indexes2]) # Flat-top in the distribution!
        # print('stdpeakIPDF = ', stdpeakIPDF)
        # print 'Peak Fraction = ', 100*peakICDF[indexes][-1]-peakICDF[indexes][0], stdpeakIPDF
        beam.bin_time()
        t = 1e12*(beam.t - np.mean(beam.t))
        dt = 0.05*(max(t) - min(t))

        ax[0][0].clear()
        ax[0][1].clear()
        ax[1][0].clear()
        ax[1][1].clear()
        ax[2][0].clear()
        ax[2][1].clear()

        exponent = np.floor(np.log10(np.abs(beam.slice_length)))
        x = 10**(12) * np.array((beam.slice_bins - np.mean(beam.t)))
        ax[0][0].plot(x, beam.slice_peak_current)
        ax[0][0].set(ylabel='I [A]')

        ax[1][0].set_xlim(min(t) - dt, max(t) + dt)
        t = 1e12*(beam.t - np.mean(beam.t))
        p = 1e-6*beam.cp
        ymax = max(p)+1
        ymin = min(p)-1
        if ymax - ymin < 5:
            ymax = np.mean(p) + 2.5
            ymin = np.mean(p) - 2.5
        ax[1][0].hist2d(t, p, bins=(50,50), cmap=plt.cm.jet, range=[[min(t), max(t)],[ymin, ymax]])
        # ax[0][2].set_ylim(top=ymax, bottom=ymin)
        ax[1][0].set(ylabel='cp [Mev]')

        ax[2][0].plot(x, 1e6*beam.slice_normalized_horizontal_emittance)
        ax[2][0].plot(x, 1e6*beam.slice_normalized_vertical_emittance)
        ax[2][0].set_ylim(top=3, bottom=0)
        ax[2][0].set(ylabel='emit_x / emit_y [m]', xlabel='t [ps]')
        ax[2][0].grid()

        ax[0][1].plot(twiss.elegant['s'], 1e12*twiss['sigma_z']/2.99e8)
        ax[0][1].set(ylabel='sigma_l [ps]')
        ax[0][1].set_ylim(bottom=0)
        print('bunch length = ', 1e12*twiss['sigma_z'][-1]/2.99e8)
        print('chirp = ', 1e-6*(max(beam.cp) - min(beam.cp)))
        # ax[0][1].fill_between(t_grid[indexes], peakICDF[indexes], 0, facecolor='gray', edgecolor='gray', alpha=0.4)

        ax[1][1].plot(twiss.elegant['s'], 0.511*twiss.elegant['pCentral0'])
        # ax[1][1].set_ylim(top=1100, bottom=0)
        ax[1][1].set(ylabel='Momentum [MeV/c]')
        ax[1][1].grid()

        ax[2][1].plot(twiss.elegant['s'], 1e3*twiss['sigma_x'])
        ax[2][1].plot(twiss.elegant['s'], 1e3*twiss['sigma_y'])
        ax[2][1].set(ylabel='sigma_x / sigma_y [m]', xlabel='s [m]')
        ax[2][1].grid()

    i = 0
    j = 0
    while i < 2:
        j = 0
        while j < 2:
            plt.setp(ax[i][j].get_xticklabels(), visible=False)
            j = j + 1
        i = i + 1
    # plt.setp(ax_yDist.get_yticklabels(), visible=False)

    plt.show()


if __name__ == '__main__':
   main()
