import os, sys
import numpy as np
import random
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.ClassFiles.genesisBeamFileMOGA import MOGA
import multiprocessing

nGenerations = 200
populationSize = 24
nChildren = 48
crossoverprobability = 0.6
mutationprobability = 0.2
ngenerations = 200

class Short1GeVMOGA(MOGA):

    parameter_names = [
        # ['CLA-L01-CAV', 'field_amplitude'],
        # ['CLA-L01-CAV', 'phase'],
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
        super(Short1GeVMOGA, self).__init__()
        self.CLARA_dir = os.path.relpath(__file__+'/../../../../CLARA/')
        self.scaling = 6
        self.sample_interval=2**(3*1)
        self.base_files = '../../../../CLARA/basefiles_6/'
        self.genesis_file = 'xara_4nm_td.in'
        self.verbose = False
        self.alphax = -0.08050797
        self.betax = 4.58986
        self.alphay = -0.05371444
        self.betay = 2.52698

def optfunc(inputargs, **kwargs):
    return moga.MOGAoptFunc(inputargs, **kwargs)

moga = Short1GeVMOGA()
moga.set_changes_file(['./nelder_mead_best_changes.yaml', '../transverse_best_changes.yaml'])
moga.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
moga.base_files = '../../../../CLARA/basefiles_6/'
moga.start_lattice = 'CLARAX'
moga.runElegant = True
moga.runGenesis = True

if __name__ == "__main__":
    if not os.name == 'nt':
        nProc = 1
    else:
        nProc = 1
    moga.create_toolbox()
    moga.create_fitness_function(optfunc)
    moga.create_weights_function(weights=(-1.0, 1.0, -1.0, 1.0, ))
    moga.create_uniform_mating_function(probability=0.3)
    moga.create_gaussian_mutation_function(probability=0.3, mu=0, sigma=[5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5,  5e6, 5,  5e6, 5, 0.005, 0.1])
    # moga.add_bounds(MIN, MAX)
    moga.create_NSGA2_selection_function()
    best = moga.load_best('./nelder_mead_best_changes.yaml')
    moga.initialise_population(best, populationSize)
    if nProc > 1:
        pool = multiprocessing.Pool(processes=nProc)
        moga.toolbox.register("map", pool.map)
    moga.initialise_MOGA(seed=5685)
    moga.eaMuPlusLambda(populationSize, nChildren, crossoverprobability, mutationprobability, ngenerations)
