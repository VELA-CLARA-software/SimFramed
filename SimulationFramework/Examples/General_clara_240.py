import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()

scaling = 3

npart=2**(3*scaling)

framework = Framework('C2V', overwrite=True)

''' You need to define the location of your astra command!'''
''' This is the default for windows...'''
# framework.astra.defineASTRACommand(['../../MasterLattice/ASTRA/astra.exe'])

framework.loadSettings('Lattices/CLA400-C2V-INJ.def')
print framework.elementIndex('CLA-L01-CAV-SOL-01')
print framework.previousElement('CLA-L01-CAV-SOL-01')
print framework.nextElement('CLA-L01-CAV-SOL-01')
# framework.astra.createInitialDistribution(npart=npart,charge=250)
# framework.createInputFiles()
# framework.runInputFiles()
beam.read_astra_beam_file('C2V/injector400.0337.001')
print 'ref_particle = ', beam.beam['reference_particle']
beam.write_HDF5_beam_file('C2V/injector400.0337.hdf5')
