import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.Framework import *
from SimulationFramework.Modules.constraints import *
import SimulationFramework.Modules.read_beam_file as rbf
import SimulationFramework.Modules.read_twiss_file as rtf
import shutil
import uuid
from functools import partial

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

class TemporaryDirectory(object):
    """Context manager for tempfile.mkdtemp() so it's usable with "with" statement."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def tempname(self):
        return 'tmp'+str(uuid.uuid4())

    def __enter__(self, dir=os.getcwd()):
        exists = True
        while exists:
            self.name = dir + '/' + self.tempname()
            if not os.path.exists(self.name):
                exists=False
                os.makedirs(self.name)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        shutil.rmtree(self.name)

class fitnessFunc():

    def __init__(self, args, tempdir, scaling=4, overwrite=True, verbose=False, summary=False, post_injector=True, parameterfile=[]):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.scaling = scaling
        self.tmpdir = tempdir
        self.verbose = verbose
        self.summary = summary
        self.overwrite = overwrite
        self.post_injector = post_injector
        self.parameters = list(args)
        self.dirname = os.path.basename(self.tmpdir)
        self.framework = Framework(self.dirname, clean=False)
        self.framework.loadSettings('Lattices/clara400_v12_v3.def')
        if isinstance(parameterfile, (list, tuple)):
            for p in parameterfile:
                self.framework.loadParametersFile(p)
        elif parameterfile is not None:
            self.framework.loadParametersFile(parameterfile)
        if not os.name == 'nt':
            self.framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
            self.framework.generator.number_of_particles = 2**(3*scaling)
        else:
            self.framework.generator.number_of_particles = 2**(3*3)
        self.framework.setElementType('quadrupole','k1', self.parameters)

    def between(self, value, minvalue, maxvalue, absolute=True):
        if absolute:
            result = max([minvalue,min([maxvalue,abs(value)])])
        else:
            result = np.sign(value)*max([minvalue,min([maxvalue,abs(value)])])
        return result

    def calculateBeamParameters(self):

        try:
            twiss = self.twiss
            if self.overwrite:
                startS = self.framework['L02'].startObject['position_start'][2]
                if self.post_injector:
                    self.framework['L02'].file_block['input']['prefix'] = '../basefiles_'+str(self.scaling)+'/'
                    self.framework.track(startfile='L02', endfile='S07')#startfile='FMS')
                else:
                    self.framework.track(endfile='S07')

            constraintsList = {}
            constraintsListQuads = {
                'max_k': {'type': 'lessthan', 'value': [abs(p) for p in self.parameters], 'limit': 2.5, 'weight': 10},

            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)
            # twiss.read_astra_emit_files( [ self.dirname+'/'+n+'.Zemit.001' for n in self.framework.fileSettings.keys() if self.framework.fileSettings[n]['code'].upper() == 'ASTRA'] )
            # constraintsListSigmas = {
            #     'max_xrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 10},
            #     'max_yrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 10},
            #     'min_xrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.1, 'weight': 10},
            #     'min_yrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.1, 'weight': 10},
            #     'last_exn': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.6, 'weight': 1},
            #     'last_eyn': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.6, 'weight': 1},
            # }
            # constraintsList = merge_two_dicts(constraintsList, constraintsListSigmas)
            twiss.read_astra_emit_files(self.dirname+'/S07.Zemit.001')
            tdc_position = self.framework['CLA-S07-TDC-01-R']['position_start'][2]
            tdc_screen_position = self.framework['CLA-S07-DIA-SCR-03-W']['position_start'][2]
            dechirper_position = self.framework['CLA-S07-DCP-01']['position_start'][2]
            constraintsListS07 = {
                'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(tdc_screen_position,'muy') - twiss.interpolate(tdc_position,'muy'), 'limit': 0.25, 'weight': 1},
                'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', tdc_position, tdc_screen_position), 'limit': 5, 'weight': 1},
                'dechirper_sigma_x': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_y': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_xy': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y') - 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.0, 'weight': 20},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListS07)
            fitness = self.cons.constraints(constraintsList)
            if self.verbose:
                print self.cons.constraintsList(constraintsList)
            # if self.summary:
            #     self.astra.createHDF5Summary(reference='Transverse_GA')
            return fitness
        except:
            return 1e6
