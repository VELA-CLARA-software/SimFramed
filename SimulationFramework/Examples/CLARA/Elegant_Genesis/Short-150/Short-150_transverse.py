import os, sys
sys.path.append('../../../../../')
import SimulationFramework.Framework as fw
# import numpy as np
from SimulationFramework.Modules.nelder_mead import nelder_mead
sys.path.append('./../../')
from Elegant.Optimise_transverse import Optimise_transverse
from ruamel import yaml

framework = fw.Framework(None)
framework.loadSettings('Lattices/clara400_v12_v3.def')
parameters = framework.getElementType('quadrupole','k1l')
names = framework.getElementType('quadrupole','objectname')
best = parameters
scaling=4

best = [-0.22789817,  0.04384427,  0.07237042,  0.11319594,  0.11633546, -0.10393746,
  0.09247306,  0.03135896, -0.06080841, -0.04500804,  0.02695322,  0.0206167,
  0.03058594, -0.04103264, -0.12178037,  0.30441347,  0.22639786, -0.08582796,
  0.37171019, -0.13076231, -0.33075536, -0.10981188,  0.43603262,  0.01990002,
  0.27849027,  0.37018414, -0.1794828,   0.03491095, -0.17280846]

with open('transverse_best_changes.yaml', 'r') as infile:
    data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
    best = [data[n]['k1l'] for n in names[:-2]]


if __name__ == "__main__":
        fit = Optimise_transverse()
        fit.setChangesFile('./nelder_mead_best_changes.yaml')
        fit.verbose = False
        fit.Nelder_Mead(best, step=0.5)
