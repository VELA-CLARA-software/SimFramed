from Framework import *
import numpy as np
import read_beam_file as rbf

beam = rbf.beam()

scaling = 4

npart=2**(3*scaling)

grids = getGrids(npart=npart)
# print grids.getGridSizes()


framework = Framework('2', overwrite=True)

framework.loadSettings('clara400_v12.def')
framework.astra.createInitialDistribution(npart=npart,charge=250)

gun = framework.getElement('CLA-HRG1-GUN-CAV')
# gun['phase'] = -5
linac = framework.getElement('CLA-L01-CAV')
# linac['phase'] = 0
sol = framework.getElement('CLA-LRG1-GUN-SOL')
# sol['field_amplitude'] = 0.35
if not os.name == 'nt':
    framework.astra.defineASTRACommand(['mpiexec','-np',str(scaling*3),'/opt/ASTRA/astra_MPICH2.sh'])
    framework.CSRTrack.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(scaling*3),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    # gun['field_definition']= '/home/jkj62/Data_Files/HRRG_1D_RF.dat'
    # sol['field_definition']= '/home/jkj62/Data_Files/HRRG_combined_sols_100mm_onaxis.dat'
    # linac['field_definition']= '/home/jkj62/Data_Files/TWS_S-DL.dat'

beam.read_astra_beam_file('2'+'/S07.4928.001')
print beam.chirp
# framework.createInputFiles()
# framework.runInputFiles()
