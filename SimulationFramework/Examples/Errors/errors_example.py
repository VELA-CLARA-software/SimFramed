import warnings
warnings.filterwarnings("ignore")
import sys, os
sys.path.append('../../../')
from SimulationFramework.Framework import *
import read_cavity_phases as rcp
from scipy.stats import truncnorm
from collections import defaultdict
import numpy as np

def create_base_file(settings='./clara400_v12_v3.def'):
    cavity_absolute_phases, cavity_absolute_momenta = rcp.get_Cavity_Phases(dir='../basefiles_4', settings='./clara400_v12_v3.def')
    lattice = Framework('test', clean=False, verbose=False)
    lattice.loadSettings(settings)
    if not os.name == 'nt':
        scaling = 4
        lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
        lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        lattice.generator.number_of_particles = 2**(3*scaling)
    else:
        lattice.generator.number_of_particles = 2**(3*3)

    lattice.defineElegantCommand(['elegant'])
    for lat in lattice.latticeObjects:
        print 'lattice = ', lat
        if hasattr(lattice[lat], 'headers') and not lat == 'generator':
            lattice.setSubDirectory('basefiles')
            lattice[lat].headers['newrun']['auto_phase'] = True
            lattice[lat].headers['newrun']['Track_All'] = False
            lattice[lat].file_block['input']['prefix'] = '../test/'
            lattice.track(files=[lat], postprocess=False)
            lattice.setSubDirectory('test')

            cavities = lattice[lat].getElementType('cavity')
            if len(cavities) > 0:
                absolute_phases, absolute_momenta = rcp.get_Cavity_Phase('basefiles/'+lat+'.log')
                print zip(cavities, absolute_phases)
                for cav, absphase in zip(cavities, absolute_phases):
                    cav['phase'] = absphase + cav['phase']
            lattice[lat].headers['newrun']['auto_phase'] = False
            lattice[lat].headers['newrun']['Track_All'] = True
        lattice.track(files=[lat], postprocess=True)

# create_base_file()

class Framework_Errors(object):

    def __init__(self, lattice=None):
        super(Framework_Errors, self).__init__()
        self.lattice = lattice
        self.errors = defaultdict(dict)

    def gaussian_error(self, loc, amplitude, cutoff=3):
        return truncnorm.rvs(a=-cutoff, b=cutoff, loc=loc, scale=amplitude)

    def create_error(self, centre, amplitude, fractional=False, type='gaussian', cutoff=3):
        if hasattr(self, type+'_error'):
            loc = 0
            error = getattr(self, type+'_error')(loc, amplitude, cutoff)
        else:
            error = 0
        return centre * error if fractional else error

    def add_error(self, element, parameter, amplitude, type='gaussian', cutoff=3, fractional=False):
        error = self.create_error(getattr(self.lattice[element], parameter), amplitude, type=type, cutoff=cutoff, fractional=fractional)
        self.lattice[element]['d'+parameter] = error
        self.errors[element]['d'+parameter] = error

    def add_position_error(self, element, dx=0, dy=0, dz=0, type='gaussian', cutoff=3, fractional=False):
        self.add_error(element, 'x', dx, type=type, cutoff=cutoff, fractional=fractional)
        self.add_error(element, 'y', dy, type=type, cutoff=cutoff, fractional=fractional)
        self.add_error(element, 'z', dz, type=type, cutoff=cutoff, fractional=fractional)

# cavity_absolute_phases, cavity_absolute_momenta =  rcp.get_Cavity_Phases(dir='./basefiles', settings='./clara400_v12_v3.def')

def load_lattice():
    global lattice
    lattice = Framework('example', clean=False, verbose=False)
    lattice.loadSettings('split.def')
    if not os.name == 'nt':
        scaling = 4
        lattice.defineASTRACommand(['mpiexec','-np',str(4*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
        lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        lattice.generator.number_of_particles = 2**(3*scaling)
    else:
        lattice.generator.number_of_particles = 2**(3*3)
    lattice.defineElegantCommand(['elegant'])

def track_with_absolute_phases(global_errors=False):
    for lat in lattice.latticeObjects:
        if hasattr(lattice[lat], 'headers') and not lat == 'generator':
            cavities = lattice[lat].getElementType('cavity')
            if len(cavities) > 0:
                absolute_phases, absolute_momenta = rcp.get_Cavity_Phase('basefiles/'+lat+'.log')
                for cav, absphase in zip(cavities, absolute_phases):
                    cav['phase'] = absphase + cav.phase
            lattice[lat].headers['newrun']['auto_phase'] = False
            lattice[lat].headers['newrun']['Track_All'] = True
            lattice[lat].headers['global_errors']['global_errors'] = global_errors
        lattice.track(files=[lat])#, preprocess=True, track=False, postprocess=False)

def vary_charge(q):
    load_lattice()
    lattice.setSubDirectory('charge_'+str(np.round(q, decimals=2)))
    lattice['generator'].charge = 250e-12 * q
    print 'Setting charge to ', lattice['generator'].charge
    track_with_absolute_phases()

def vary_initial_t(t):
    load_lattice()
    lattice.setSubDirectory('initial_t_'+str(np.round(t, decimals=2))+'ps')
    lattice['injector400'].headers['newrun'].toff = 1e-3 * t
    print 'Setting T offset to ', 1e3*lattice['injector400'].headers['newrun'].toff, 'ps'
    track_with_absolute_phases()

def vary_initial_x(x):
    load_lattice()
    lattice.setSubDirectory('initial_x_'+str(np.round(x, decimals=2))+'mm')
    lattice['injector400'].headers['newrun'].xoff = x
    print 'Setting X offset to ', lattice['injector400'].headers['newrun'].xoff, 'mm'
    track_with_absolute_phases()

def vary_initial_y(y):
    load_lattice()
    lattice.setSubDirectory('initial_y_'+str(np.round(y, decimals=2))+'mm')
    lattice['injector400'].headers['newrun'].yoff = y
    print 'Setting Y offset to ', lattice['injector400'].headers['newrun'].yoff, 'mm'
    track_with_absolute_phases()

def vary_quad_dx(quad, x):
    load_lattice()
    if lattice.getElement(quad) is not None:
        lattice.setSubDirectory(quad+'_x_'+str(np.round(x, decimals=2))+'mm')
        lattice[quad].dx = 1e-3*x
        print 'Setting ', quad, 'X offset to ', 1e3*lattice[quad].dx, 'mm'
        track_with_absolute_phases()

def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step

# for q in frange(0.91,1.1,0.01):
#     vary_charge(q)

# toff is in ps!
# for t in frange(-1,1,0.1):
#     vary_initial_t(t)

# xoff is in mm!
# for x in frange(-1,1,0.1):
#     vary_initial_x(x)
# for y in frange(-1,1,0.1):
#     vary_initial_y(y)
load_lattice()
quads = [a.objectName for a in lattice.getElementType('Quadrupole')]
# print quads
for q in quads[1:]:
    for x in frange(-1,1.01,1):
        vary_quad_dx(q, x)
# errors = Framework_Errors(lattice)
# print lattice['CLA-S02-MAG-QUAD-01'].start
# errors.add_position_error('CLA-S02-MAG-QUAD-01', dx=1e-3, dy=1e-3, dz=1e-3, fractional=False)
exit()
