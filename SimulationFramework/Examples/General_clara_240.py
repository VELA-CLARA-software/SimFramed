import sys, os, math
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()
degree = 1./360.*2.*math.pi

scaling = 4

npart=2**(3*scaling)


''' ############# CLARA 400 to VELA #############'''

# framework = Framework('C2V', overwrite=True)
#
# framework.loadSettings('Lattices/CLA400-C2V-INJ.def')
# framework.astra.createInitialDistribution(npart=npart,charge=250)
# framework.createRunProcessInputFiles(files=['vela'])


''' ############# CLARA 400 #############'''

framework = Framework('CLARA', overwrite=True)

framework.loadSettings('Lattices/clara400_v12.def')
framework.astra.createInitialDistribution(npart=npart,charge=250)
framework.createRunProcessInputFiles(run=True)#, files=['VBC'])
