import os, sys
import numpy as np
import random
sys.path.append('./../')
from genesisBeamFileMOGA import MOGA
import multiprocessing

nGenerations = 200
populationSize = 24
nChildren = 48
crossoverprobability = 0.6
mutationprobability = 0.2
ngenerations = 200

# best = [31295498.77911233, -21.796417988221847, 28533996.939464197, -4.3535424459076175, 25118020.516425118,
# 186.03322780846452, 31693858.476702053, 49.0649623919314, -0.12195857313414843, 0.9046648051503599]
best = [31384565,-23.3,29675255,-3.4000000000000004,27484288,187.5,31191886,50.7,-0.12240000000000001,0.6128]
MIN = [0, -90, 0, -90, 0, 90, 0, -90, -0.2, 0.0]
MAX = [33e6, 90, 33e6, 90, 45e6, 270, 32e6, 90, -0.05, 3]

def optfunc(*args, **kwargs):
    return moga.MOGAoptFunc(*args, **kwargs)

moga = MOGA()

if __name__ == "__main__":
    if not os.name == 'nt':
        nProc = 8
    else:
        nProc = 2
    moga.create_toolbox()
    moga.create_fitness_function(optfunc, scaling=5, post_injector=True, changes='./transverse_best_changes.yaml', lattice='./Short-240.def')
    moga.create_weights_function(weights=(-1.0, 1.0, -1.0, 1.0, ))
    moga.create_uniform_mating_function(probability=0.3)
    moga.create_gaussian_mutation_function(probability=0.3, mu=0, sigma=[1e6,2,1e6,2,2e6,2,1e6,2,0.01,0.2])
    # moga.add_bounds(MIN, MAX)
    moga.create_NSGA2_selection_function()
    moga.initialise_population(best, populationSize)
    if nProc > 1:
        pool = multiprocessing.Pool(processes=nProc)
        moga.toolbox.register("map", pool.map)
    moga.initialise_MOGA(seed=5685)
    moga.eaMuPlusLambda(populationSize, nChildren, crossoverprobability, mutationprobability, ngenerations)
