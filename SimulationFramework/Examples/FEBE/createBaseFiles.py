import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../'))
from SimulationFramework.Framework import *

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling, changes=None):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('FEBE.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(4*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(4*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.generator.number_of_particles = 2**(3*scaling)
    if not changes is None:
        framework.load_changes_file(changes)
    framework.track()

## Modify as appropriate! ##
for i in [5]:
    create_base_files(i, './nelder_mead_best_changes.yaml')
exit()
