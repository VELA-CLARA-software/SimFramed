import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
import SimulationFramework.Modules.read_twiss_file as rtf
import numpy as np
ncpu = 10
beam = rbf.beam()
twiss = rtf.twiss()
import csv
################################  CSRTrack #####################################

# framework = Framework('VBC_CSRTrack')
# if not os.name == 'nt':
#     framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
#     framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
#     framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
# framework.defineElegantCommand(['elegant'])
#
# framework.loadSettings('Lattices/clara400_v12_v3.def')
# framework['VBC'].file_block['input']['prefix'] = '../basefiles_5/'
# framework.track(startfile='VBC', endfile='S07')


################################  ELEGANT ######################################

# framework = Framework('VBC_Elegant')
# if not os.name == 'nt':
#     framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
#     framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
#     framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
# framework.defineElegantCommand(['elegant'])
#
# framework.loadSettings('Lattices/clara400_v12_v3.def')
# framework.change_Lattice_Code('VBC', 'elegant')
# framework['VBC'].file_block['input']['prefix'] = '../basefiles_5/'
# framework.track(startfile='VBC', endfile='S07')

################################  ELEGANT ######################################

# framework = Framework('CLA10BA1', clean=True)
# if not os.name == 'nt':
#     framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
#     framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
#     framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
# framework.defineElegantCommand(['elegant'])
#
# framework.loadSettings('Lattices/CLA10-BA1.def')
# framework.change_Lattice_Code('All','elegant')
# framework['S02'].prefix = '../../basefiles_4/'
# framework.track(startfile='S02', endfile='BA1_dipole')

basedir = 'dipole_scan'
framework = Framework(basedir, clean=True, verbose=False)
if not os.name == 'nt':
    framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
    framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
    framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
framework.defineElegantCommand(['elegant'])

framework.loadSettings('Lattices/CLA10-BA1.def')
framework.change_Lattice_Code('All','elegant')
framework['S02'].prefix = '../../basefiles_4/'
etax = []
etax2 = []
kValues = np.arange(0,3,0.1)
for k in kValues:#np.arange(0,1,0.1):
    print('setting k = ', k)
    framework.modifyElement('EBT-BA1-MAG-QUAD-07','k1l', k)
    framework.track(startfile='S02', endfile='BA1_dipole')
    twiss.reset_dicts()
    twiss.read_sdds_file( basedir+'/BA1_dipole.mat' )
    etax.append(twiss.elegant['R16'][-1])
    etax2.append(twiss.elegant['T166'][-1])
    print(k, twiss.elegant['R16'][-1], twiss.elegant['T166'][-1], twiss.elegant['R16'][-1] / twiss.elegant['T166'][-1])

data = zip(kValues, etax, etax2, np.array(etax)/np.array(etax2))
with open('dispersion.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(data)
