import os, sys
from copy import copy
import numpy as np
import random
sys.path.append('./../../../')
from elegantGA import GA
from XARA_nelder_mead import XARA
import multiprocessing
from SimulationFramework.Examples.CLARA.Elegant.Optimise_longitudinal_Elegant import saveState
import SimulationFramework.Modules.id_number as idn
import SimulationFramework.Modules.id_number_server as idnserver
import best_value as bestclient
import best_value_server as bestserver
from shutil import copyfile

nGenerations = 200
populationSize = 24
nChildren = 48
crossoverprobability = 0.6
mutationprobability = 0.2
ngenerations = 200

class XARA_GA(XARA):
    def __init__(self):
        super(XARA_GA, self).__init__()
        self.subdir = 'GA'
        self.optdir = self.subdir + '/iteration_'
        self.best_changes = './GA_best_changes.yaml'
        self.bestfit = 1e12

    def OptimisingFunction(self, inputargs, *args, **kwargs):
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
        self.opt_iteration += 1
        constraintsList = self.calculate_constraints()
        fitvalue = self.cons.constraints(constraintsList)
        print('fitvalue[', self.opt_iteration-1, '] = ', fitvalue)
        saveState(self.subdir, inputargs, self.opt_iteration-1, fitvalue)
        bestcli.set_best(fitvalue, self.opt_iteration)
        self.bestfit = bestcli.get_best()
        # print('best = ', best)
        if fitvalue <= self.bestfit:
            print(self.cons.constraintsList(constraintsList))
            print('!!!!!!  New best = ', fitvalue, inputargs)
            copyfile(dir+'/changes.yaml', self.best_changes)
            self.bestfit = fitvalue
        return (fitvalue, )

def optfunc(*args, **kwargs):
    return opt.optfunc(*args, **kwargs)

ga = GA()
opt = XARA_GA()
opt.set_changes_file(['nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
opt.set_start_file('CLARAX')
best = opt.load_best('nelder_mead_best_changes.yaml')

if __name__ == "__main__":
    server = idnserver.zmqServer()
    server.daemon = True
    server.start()

    bestsrv = bestserver.zmqServer()
    bestsrv.daemon = True
    bestsrv.start()

    if not os.name == 'nt':
        nProc = 8
    else:
        nProc = 1
        opt.plotting = True
    ga.create_toolbox()
    ga.create_fitness_function(optfunc)
    ga.create_weights_function(weights=(-1.0, ))
    ga.create_uniform_mating_function(probability=0.3)
    ga.create_gaussian_mutation_function(probability=0.3, mu=0, sigma=[5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5,  5e6, 5,  5e6, 5, 0.005, 0.1])
    # ga.add_bounds(MIN, MAX)
    ga.create_tournament_selection_function(tournsize=3)
    ga.initialise_population(best, populationSize)
    if nProc > 1:
        pool = multiprocessing.Pool(processes=nProc)
        ga.toolbox.register("map", pool.map)
    ga.initialise_GA(seed=5685)
    ga.eaMuPlusLambda(populationSize, nChildren, crossoverprobability, mutationprobability, ngenerations)
