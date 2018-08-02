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

lattice = Framework('short_240', clean=False, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
if not os.name == 'nt':
    scaling = 5
    lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
    # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
    lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    lattice.defineElegantCommand(['elegant'])
    lattice.generator.number_of_particles = 2**(3*scaling)
else:
    lattice.generator.number_of_particles = 2**(3*3)

lattice['CLA-HRG1-GUN-SOL']['field_amplitude'] = 0.345
lattice['CLA-L01-CAV']['field_amplitude'] = 1e6*27.26
lattice['CLA-L01-CAV']['phase'] = -20.13
lattice['CLA-L02-CAV']['field_amplitude'] = 1e6*23.685
lattice['CLA-L02-CAV']['phase'] = -18.342
lattice['CLA-L03-CAV']['field_amplitude'] = 1e6*18.5463
lattice['CLA-L03-CAV']['phase'] = -10.319
lattice['CLA-L4H-CAV']['field_amplitude'] = 1e6*22.91
lattice['CLA-L4H-CAV']['phase'] = 181.5
lattice['CLA-L04-CAV']['field_amplitude'] = 1e6*31.522
lattice['CLA-L04-CAV']['phase'] = 6.032
lattice['bunch_compressor'].update(dipoleangle=0.11715)
lattice['CLA-L01-WAKE'].scale_kick = 0
lattice['CLA-L02-WAKE'].scale_kick = 0
lattice['CLA-L03-WAKE'].scale_kick = 0
lattice['CLA-L04-WAKE'].scale_kick = 0
lattice.track(track=True, postprocess=True)
