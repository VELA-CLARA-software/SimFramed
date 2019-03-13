import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *

# lattice = framework('C2V', clean=False, verbose=False)
# lattice.loadSettings('Lattices/cla400-ba1.def')
# if not os.name == 'nt':
#     scaling = 5
#     lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
#     lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
#     lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
#     lattice.generator.number_of_particles = 2**(3*scaling)
# else:
#     lattice.generator.number_of_particles = 2**(3*3)
# lattice.track()#startfile='C2V')

lattice = Framework('example', clean=False, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
if not os.name == 'nt':
    scaling = 5
    lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
    # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
    lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    lattice.generator.number_of_particles = 2**(3*scaling)
else:
    lattice.generator.number_of_particles = 2**(3*3)
lattice.defineElegantCommand(['elegant'])
print lattice['S02'].getElement('CLA-S02-MAG-QUAD-01','k1')
lattice.modifyElement('CLA-S02-MAG-QUAD-01','k1', 10)
print lattice['S02'].getElement('CLA-S02-MAG-QUAD-01','k1')
lattice.change_Lattice_Code('S02','elegant')
lattice.track(files=['S02'])#,preprocess=True, track=False, postprocess=False)
