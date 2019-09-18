import os, sys
import numpy as np
import random
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.ClassFiles.optimise_longitudinal_Elegant_MOGA import MOGA
import multiprocessing
import traceback
import SimulationFramework.Modules.id_number as idn
import SimulationFramework.Modules.id_number_server as idnserver
from copy import copy
import deap, deap.base, deap.creator

populationSize = 16
nChildren = 2*populationSize
crossoverprobability = 0.6
mutationprobability = 0.2
ngenerations = 10

class UltrashortMOGA(MOGA):

    parameter_names = [
        ['CLA-HRG1-GUN-CAV', 'phase'],
        ['CLA-HRG1-GUN-SOL', 'field_amplitude'],
        ['CLA-L01-CAV-SOL-01', 'field_amplitude'],
        ['CLA-L01-CAV-SOL-02', 'field_amplitude'],
        ['CLA-L01-CAV', 'field_amplitude'],
        ['CLA-L01-CAV', 'phase'],
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV', 'field_amplitude'],
        ['CLA-L04-CAV', 'phase'],
        #if add/remove parameters here, change mutation function sigmas below
    ]

    def __init__(self):
        super(UltrashortMOGA, self).__init__()
        # self.CLARA_dir = os.path.relpath(__file__+'/../../../../CLARA/')
        self.scaling = 3
        #self.sample_interval=2**(3*1)
        # self.base_files = '../../../../CLARA/basefiles_6/'
        self.verbose = False
        self.ncpu = 1
        self.change_to_elegant = False
        self.post_injector = False
        self.doTracking = True

    def before_tracking(self):
        elements = self.framework.elementObjects.values()
        for e in elements:
            e.lsc_enable = True
            e.lsc_bins = 200
            e.current_bins = 0
            e.longitudinal_wakefield_enable = True
            e.transverse_wakefield_enable = True
        lattices = self.framework.latticeObjects.values()
        for l in lattices:
            l.lscDrifts = True
            l.lsc_bins = 200
        #add TURN VBC OFF HERE
        #set QUADS HERE

    def MOGAoptFunc(self, inputargs, *args, **kwargs):
        energy_spread, peak_current = self.OptimisingFunction(inputargs, **kwargs)
        fitness = energy_spread/peak_current #doesn't do much atm
        print ('fitvalue[', self.opt_iteration, '] - eSpread=', energy_spread, '  peakI=', peak_current, '  fitness=', fitness)
        self.saveState(inputargs, self.opt_iteration, [energy_spread, peak_current], fitness)
        return energy_spread, peak_current

    def OptimisingFunction(self, inputargs, **kwargs):
        self.optdir = 'MOGA/iteration_'
        parameternames = copy(self.parameter_names)

        self.inputlist = list(map(lambda a: a[0]+[a[1]], zip(parameternames, inputargs)))

        self.linac_fields = np.array([i[2] for i in self.inputlist if i[1] == 'field_amplitude'])
        self.linac_phases = np.array([i[2] for i in self.inputlist if i[1] == 'phase'])

        idclient = idn.zmqClient()
        n =  idclient.get_id()
        self.opt_iteration = n

        dir = self.optdir+str(self.opt_iteration)
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.setup_lattice(self.inputlist, dir)
        self.before_tracking()
        fitvalue = self.track(endfile='S07')

        self.beam.read_HDF5_beam_file(self.dirname+'/CLA-S07-MARK-03.hdf5')
        self.beam.slices = 50
        self.beam.bin_time()
        sigmat = 1e12*np.std(self.beam.t)
        sigmap = np.std(self.beam.cp)
        meanp = np.mean(self.beam.cp)
        emitx = 1e6*self.beam.normalized_horizontal_emittance
        emity = 1e6*self.beam.normalized_horizontal_emittance
        fitp = 100*sigmap/meanp #percentage energy spread
        peakI, peakIstd, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum, peakIDensity = self.beam.sliceAnalysis()

        if isinstance(self.opt_iteration, int):
            self.opt_iteration += 1

        return fitp, peakI

def optfunc(inputargs, **kwargs):
    try:
        return moga.MOGAoptFunc(inputargs, **kwargs)
    except Exception as e:
        print(traceback.format_exc())
        return (1, 0)

moga = UltrashortMOGA()
moga.set_lattice_file('Lattices/clara400_v12_v3.def')
moga.start_lattice = 'generator'

if os.name == 'nt':
    nProc = 16 #change this, if more than 1 doesn't show useful errors
else:
    nProc = 11
moga.create_toolbox()
moga.create_fitness_function(optfunc)
moga.create_weights_function(weights=(-1.0, 1.0, )) #minus is minimise
moga.create_uniform_mating_function(probability=0.3)
# gun phase, gun sol, l01 sol1, l01 sol2, l01 amp, l01 phase, l02 amp, l02 phase, l03 amp, l03 phase, l4h amp, l4h phase, l04 amp, l04 phase
moga.create_gaussian_mutation_function(probability=0.3, mu=0, sigma=[1, 0.01, 0.01, 0.01, 1e6, 1, 1e6, 1, 1e6, 1, 1e6, 1, 1e6, 1])
# moga.add_bounds(MIN, MAX)
moga.create_NSGA2_selection_function()

if __name__ == "__main__":
    best = moga.load_best('./GA_best_changes.yaml')
    moga.initialise_population(best, populationSize, sigma=[5, 0.05, 0.05, 0.05, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5])

    if nProc > 1:
        pool = multiprocessing.Pool(processes=nProc)
        moga.toolbox.register("map", pool.map)
    moga.initialise_MOGA(seed=23158)
    moga.eaMuPlusLambda(populationSize, nChildren, crossoverprobability, mutationprobability, ngenerations)
