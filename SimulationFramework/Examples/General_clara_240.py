from Framework import *
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
