import sys, os
from copy import copy
import numpy as np
sys.path.append('./../../../')
from SimulationFramework.Examples.CLARA.Elegant.Optimise_longitudinal_Elegant import Optimise_Elegant, saveState
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from functools import partial
from ruamel import yaml
from shutil import copyfile
import SimulationFramework.Framework as fw
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import operator
import random
import multiprocessing
from scoop import futures
import csv
from SimulationFramework.Modules.optimiser import optimiser
eaopt = optimiser()
# sys.path.append('./../CLARA/Elegant_Genesis/')
import SimulationFramework.Modules.id_number as idn
import SimulationFramework.Modules.id_number_server as idnserver
import best_value as bestclient
import best_value_server as bestserver

# import matplotlib.pyplot as plt
# plt.ion()
# plt.show()
# fig, ax = plt.subplots(2, 3, sharey=False,
#                        figsize=(16, 6))

class XARA(Optimise_Elegant):

    parameter_names = [
        ['CLA-L01-CAV', 'field_amplitude'],
        ['CLA-L01-CAV', 'phase'],
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV-01', 'field_amplitude'],
        ['CLA-L04-CAV-01', 'phase'],
        ['CLA-L04-CAV-02', 'field_amplitude'],
        ['CLA-L04-CAV-02', 'phase'],
        ['CLA-L04-CAV-03', 'field_amplitude'],
        ['CLA-L04-CAV-03', 'phase'],
        ['CLA-L04-CAV-04', 'field_amplitude'],
        ['CLA-L04-CAV-04', 'phase'],
        ['bunch_compressor', 'angle'],
        ['CLA-S07-DCP-01', 'factor'],
    ]

    def __init__(self):
        super(XARA, self).__init__()
        self.scaling = 6
        self.fit.CLARA_dir = os.path.relpath(__file__+'/../../../../CLARA/')
        self.optfunc = partial(self.OptimisingFunction2, scaling=self.scaling, post_injector=self.POST_INJECTOR, verbose=False, sample_interval=2**(3*2))
        self.verbose = False
        self.bestfit = 1e12

    def load_best(self, filename):
        with open(filename, 'r') as infile:
            data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
            # best = [data[n]['k1l'] for n in parameter_names]
            best = []
            for n, p in self.parameter_names:
                if n in data:
                    best.append(data[n][p])
                elif n == 'bunch_compressor' and p == 'set_angle':
                    best.append(data['CLA-VBC-MAG-DIP-01']['angle'])
                else:
                    print((n, p))
                    if not hasattr(self, 'framework'):
                        self.framework = fw.Framework(None)
                        self.framework.loadSettings(self.lattice)
                    best.append(self.framework[n][p])
            self.best = best
        return best

    def OptimisingFunction2(self, inputargs, *args, **kwargs):
        try:
            idclient = idn.zmqClient()
            n =  idclient.get_id()
            self.opt_iteration = n
            bestcli = bestclient.zmqClient()
            if not self.POST_INJECTOR:
                parameternames = self.injector_parameter_names + self.parameter_names
            else:
                parameternames = copy(self.parameter_names)
            self.inputlist = [a[0]+[a[1]] for a in zip(parameternames, inputargs)]

            self.linac_fields = np.array([i[2] for i in self.inputlist if i[1] == 'field_amplitude'])
            self.linac_phases = np.array([i[2] for i in self.inputlist if i[1] == 'phase'])

            if 'dir' in list(kwargs.keys()):
                dir = kwargs['dir']
                del kwargs['dir']
            else:
                dir = self.optdir+str(self.opt_iteration)

            self.fit.setup_lattice(self.inputlist, dir, changes=self.changes, *args, **kwargs)
            fitvalue = self.fit.calculateBeamParameters()
            constraintsList = self.calculate_constraints()
            fitvalue = self.cons.constraints(constraintsList)
            print('fitvalue[', self.opt_iteration, '] = ', fitvalue)
            # saveState(self.subdir, inputargs, self.opt_iteration, fitvalue)
            # print('bestcli = ', bestcli.set_best(fitvalue, self.opt_iteration))
            if fitvalue < bestcli.set_best(fitvalue, self.opt_iteration):
                print(self.cons.constraintsList(constraintsList))
                print('!!!!!!  New best = ', fitvalue, inputargs)
                copyfile(dir+'/changes.yaml', self.best_changes)
                self.bestfit = fitvalue
        except:
            fitvalue = 1e12
        return (fitvalue,)

    def calculate_constraints(self):
        self.twiss = self.fit.twiss
        twiss = self.twiss
        self.beam = self.fit.beam
        constraintsList = {}
        quadkls = self.fit.framework.getElementType('quadrupole','k1l')
        quadlengths = self.fit.framework.getElementType('quadrupole','length')
        quadnames = self.fit.framework.getElementType('quadrupole','objectname')

        self.beam.read_HDF5_beam_file(self.fit.dirname+'/CLA-S07-APER-01.hdf5')
        self.beam.slice_length = 0.05e-12

        t = 1e12*(self.beam.t-np.mean(self.beam.t))
        t_grid = np.linspace(min(t), max(t), 2**8)
        peakIPDF = self.beam.PDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**4))
        peakICDF = self.beam.CDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**4))
        peakIFWHM, indexes = self.beam.FWHM(t_grid, peakIPDF, frac=0.1)
        peakIFWHM2, indexes2 = self.beam.FWHM(t_grid, peakIPDF, frac=20)
        stdpeakIPDF = (max(peakIPDF[indexes2]) - min(peakIPDF[indexes2]))/np.mean(peakIPDF[indexes2]) # Flat-top in the distribution!
        # print('stdpeakIPDF = ', stdpeakIPDF)
        # print 'Peak Fraction = ', 100*peakICDF[indexes][-1]-peakICDF[indexes][0], stdpeakIPDF
        self.beam.bin_time()
        t = 1e12*(self.beam.t - np.mean(self.beam.t))
        dt = 0.05*(max(t) - min(t))

        sigmat = np.std(t)
        sigmap = np.std(self.beam.p)
        meanp = np.mean(self.beam.p)
        fitp = 100*sigmap/meanp
        peakI, peakIstd, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum, peakIDensity = self.beam.sliceAnalysis()
        chirp = self.beam.chirp

        slinac_fields = np.array([1e-6*self.linac_fields[i] for i in [0,1,2]])
        x4hlinac_fields = np.array([1e-6*self.linac_fields[i] for i in [3]])
        xlinac_fields = np.array([1e-6*self.linac_fields[i] for i in [4,5,6,7]])

        twiss.read_elegant_twiss_files( self.fit.dirname+'/CLARAX.twi' )
        ipindex = list(twiss.elegant['ElementName']).index('CLA-S07-APER-01')
        constraintsListXARA = {
            # 'ip_enx': {'type': 'lessthan', 'value': 1e6*twiss.elegant['enx'][ipindex], 'limit': 2, 'weight': 0},
            # 'ip_eny': {'type': 'lessthan', 'value': 1e6*twiss.elegant['eny'][ipindex], 'limit': 0.5, 'weight': 2.5},
            'field_max_s': {'type': 'lessthan', 'value': slinac_fields, 'limit': 32, 'weight': 300},
            'field_max_x': {'type': 'lessthan', 'value': xlinac_fields, 'limit': 80, 'weight': 300},
            'field_max_x4h': {'type': 'lessthan', 'value': x4hlinac_fields, 'limit': 35, 'weight': 300},
            'momentum_max': {'type': 'lessthan', 'value': 0.511*twiss.elegant['pCentral0'][ipindex], 'limit': 1050, 'weight': 250},
            'momentum_min': {'type': 'greaterthan', 'value': 0.511*twiss.elegant['pCentral0'][ipindex], 'limit': 990, 'weight': 150},
            'peakI_min': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 950, 'weight': 200},
            'peakI_max': {'type': 'lessthan', 'value': abs(peakI), 'limit': 1050, 'weight': 200},
            # 'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.1, 'weight': 2},
            'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 0.75, 'weight': 15},
            'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 0.75, 'weight': 1.5},
            'peakIFWHM': {'type': 'lessthan','value': peakIFWHM, 'limit': 0.5, 'weight': 100},
            'stdpeakIFWHM': {'type': 'lessthan','value': stdpeakIPDF, 'limit': 1.5, 'weight': 50},
            'peakIFraction': {'type': 'greaterthan','value': 100*peakICDF[indexes][-1]-peakICDF[indexes][0], 'limit': 90, 'weight': 20},
            'chirp': {'type': 'lessthan', 'value': abs(chirp), 'limit': 1, 'weight': 25},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListXARA)

        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print((self.cons.constraintsList(constraintsList)))
        return constraintsList

opt = XARA()
opt.set_changes_file(['nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
opt.set_start_file('CLARAX')
best = opt.load_best('GA_best_changes.yaml')

opt.subdir = 'longitudinal_GA'
opt.optdir = opt.subdir + '/iteration_'
opt.best_changes = './GA_best_changes.yaml'

# startranges = [[0.95*i, 1.05*i] if abs(i) > 0 else [-20,20] for i in best]
startranges = [[b-s, b+s] for b,s in zip(best,[1e6, 1, 1e6, 1, 1e6, 1, 1e6, 1, 1e6, 1, 1e6, 1,  1e6, 1,  1e6, 1, 0.025, 0.05])]
# print(('startranges = ', startranges))
generateHasBeenCalled = False
def generate():
    global generateHasBeenCalled
    if not generateHasBeenCalled:
        generateHasBeenCalled = True
        return creator.Individual(list(best))
    else:
        return creator.Individual(random.uniform(a,b) for a,b in startranges)

# print generate()

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", generate)

# Structure initializers
toolbox.register("Individual", generate)
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

if os.name == 'nt':
    toolbox.register("evaluate", opt.optfunc)
else:
    toolbox.register("evaluate", opt.optfunc)

toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":

    server = idnserver.zmqServer()
    server.daemon = True
    server.start()

    bestsrv = bestserver.zmqServer()
    bestsrv.daemon = True
    bestsrv.start()

    global hof
    random.seed(64)

    out = open('best_solutions_longitudinal_GA_elegant.csv','w')
    hoffile='XARA_HOF_longitudinal_Elegant.csv'
    # FF.csv_out = csv.writer(out)

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=8)
    else:
        pool = multiprocessing.Pool(processes=6)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=24)
    else:
        pop = toolbox.population(n=24)

    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, logbook = eaopt.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=50,
                            stats=stats, halloffame=hof, verbose=True,
                            hoffile=hoffile)

    pool.close()
    # print 'pop = ', pop
    print(logbook)
    print(hof)

    # try:
    print('best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/XARA_best_longitudinal_elegant', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=post_inj))
