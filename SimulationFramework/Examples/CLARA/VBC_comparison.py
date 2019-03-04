import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *

ncpu = 20

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

################################  ASTRA ######################################

framework = Framework('Phase_Comparison_ASTRA')
if not os.name == 'nt':
    framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
    framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
    framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
framework.defineElegantCommand(['elegant'])

framework.loadSettings('Lattices/clara400_v12_v3.def')
# framework.change_Lattice_Code('S02', 'ASTRA')
# framework.change_Lattice_Code('L02', 'ASTRA')
# framework.change_Lattice_Code('S03', 'ASTRA')
# framework.change_Lattice_Code('L03', 'ASTRA')
# framework.change_Lattice_Code('S04', 'ASTRA')
# framework.change_Lattice_Code('L4H', 'ASTRA')
# framework.change_Lattice_Code('S05', 'ASTRA')
framework.change_Lattice_Code('VBC', 'ASTRA')
framework.change_Lattice_Code('S06', 'ASTRA')
framework.change_Lattice_Code('L04', 'ASTRA')
framework.change_Lattice_Code('S07', 'ASTRA')
framework['VBC'].file_block['input']['prefix'] = '../basefiles_5/'
framework.track(startfile='VBC', endfile='S07')


################################  ELEGANT ######################################

framework = Framework('Phase_Comparison_Elegant')
if not os.name == 'nt':
    framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
    framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
    framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
framework.defineElegantCommand(['elegant'])

framework.loadSettings('Lattices/clara400_v12_v3.def')
framework.change_Lattice_Code('VBC', 'elegant')
framework.change_Lattice_Code('S06', 'ASTRA')
framework.change_Lattice_Code('L04', 'ASTRA')
framework.change_Lattice_Code('S07', 'ASTRA')
framework['VBC'].file_block['input']['prefix'] = '../basefiles_5/'
framework.track(startfile='VBC', endfile='S07')
