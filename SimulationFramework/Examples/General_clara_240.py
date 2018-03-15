# if __name__ == '__main__':
#     if __package__ is None:
#         import sys
#         from os import path
#         print path.dirname( path.abspath(__file__) )
#         sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
#         import Framework
#     else:
#         from .. import Framework
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *
# import numpy as np
# import read_beam_file as rbf
#
# beam = rbf.beam()

scaling = 4

npart=2**(3*scaling)

grids = getGrids(npart=npart)
# print grids.getGridSizes()


framework = Framework('C2V', overwrite=True)

framework.loadSettings('../../MasterLattice/Lattices/CLA400-C2V-SP1.def')
framework.astra.createInitialDistribution(npart=npart,charge=250)
framework.createInputFiles()
# framework.runInputFiles()
