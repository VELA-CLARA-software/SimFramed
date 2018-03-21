import sys, os, math
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()
degree = 1./360.*2.*math.pi

scaling = 3

npart=2**(3*scaling)

framework = Framework('C2V', overwrite=True)

''' You need to define the location of your astra command!'''
''' This is the default for windows...'''
# framework.astra.defineASTRACommand(['../../MasterLattice/ASTRA/astra.exe'])

framework.loadSettings('Lattices/CLA400-C2V-INJ.def')
# print framework.nextElement('CLA-L01-CAV')
# print framework.elementIndex('CLA-L01-CAV')
# print framework.previousElement('CLA-L01-CAV')
framework.astra.createInitialDistribution(npart=npart,charge=250)
framework.createRunProcessInputFiles()#files=['vela'])
# framework.runInputFiles()
# framework.postProcessInputFiles()
# beam.read_HDF5_beam_file('C2V/injector400.0098.hdf5')
# print beam.beam['x'][:10]
# beam.rotate_beamXZ(180*degree)
# print beam.beam['x'][:10]
# beam.unrotate_beamXZ()
# print beam.beam['x'][:10]
# beam.read_astra_beam_file('C2V/injector400.0337.001')
# print 'ref_particle = ', beam.beam['reference_particle']
# beam.write_HDF5_beam_file('C2V/injector400.0337.hdf5')
