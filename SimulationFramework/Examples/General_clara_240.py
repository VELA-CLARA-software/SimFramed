import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *


scaling = 3

npart=2**(3*scaling)

framework = Framework('C2V', overwrite=True)

''' You need to define the location of your astra command!'''
''' This is the default for windows...'''
# framework.astra.defineASTRACommand(['../../MasterLattice/ASTRA/astra.exe'])

framework.loadSettings('Lattices/CLA400-C2V-INJ.def')
framework.astra.createInitialDistribution(npart=npart,charge=250)
framework.createInputFiles()
# framework.runInputFiles()
