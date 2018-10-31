import sys, os
sys.path.append('../../../')
from SimulationFramework.Framework import *

lattice = Framework('.', clean=False, verbose=True)
lattice.loadSettings('Lattices/CLA10-BA1.def')
if not os.name == 'nt':
    scaling = 5
    lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
    # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
    lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
lattice.defineElegantCommand(['elegant'])

# lattice['BA1'].file_block['input']['prefix'] = '../BA1/'
lattice.generator.particles = 2**(3*3)
lattice.track(startfile='S02BA1', track=True, postprocess=True, preprocess=False)
