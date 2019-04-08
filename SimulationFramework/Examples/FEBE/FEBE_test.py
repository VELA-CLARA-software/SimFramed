import sys, os
sys.path.append('../../../')
import SimulationFramework.Framework as fw
import SimulationFramework.Modules.read_twiss_file as rtf
import numpy as np

lattice = fw.Framework('example', clean=False, verbose=True)
lattice.loadSettings('FEBE.def')
if not os.name == 'nt':
    scaling = 5
    lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
    # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
    lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    lattice.generator.number_of_particles = 2**(3*scaling)
else:
    lattice.generator.number_of_particles = 2**(3*5)
lattice.defineElegantCommand(['elegant'])
lattice['S02'].prefix = '../../CLARA/basefiles_5/'
lattice.change_Lattice_Code('All','elegant')
# lattice.load_changes_file('./CLARA/Elegant_Genesis/Short-240/transverse_best_changes.yaml')
lattice.track(startfile='S02', endfile='FEBE')
