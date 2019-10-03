import numpy as np
import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
import SimulationFramework.Framework as fw
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

    def __init__(self, dir=os.getcwd(), *args, **kwargs):
        self.dir = dir
        self.args = args
        self.kwargs = kwargs

    def tempname(self):
        return 'tmp'+str(uuid.uuid4())

    def __enter__(self):
        exists = True
        while exists:
            self.name = self.dir + '/' + self.tempname()
            if not os.path.exists(self.name):
                exists=False
                os.makedirs(self.name)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.name)
        except:
            pass

class fitnessFunc(object):

    def __init__(self, lattice='Lattices/clara400_v12_v3_elegant_jkj.def'):
        self.CLARA_dir = os.path.relpath(os.path.dirname(os.path.dirname( os.path.abspath(__file__))))# if CLARA_dir is None else CLARA_dir
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.lattice_file = lattice
        self.base_files = None
        self.verbose = False
        self.start_lattice = None
        self.scaling = 6
        self.overwrite = True
        self.post_injector = True
        self.startcharge = None
        self.changes = None
        self.sample_interval = 1

    def set_CLARA_directory(self, clara_dir):
        self.CLARA_dir = clara_dir

    def set_start_file(self, file):
        self.start_lattice = file

    def setup_lattice(self, inputargs, tempdir, *args, **kwargs):
        self.dirname = tempdir
        self.input_parameters = list(inputargs)

        self.npart=2**(3*self.scaling)
        ncpu = self.scaling*3

        self.framework = fw.Framework(self.dirname, overwrite=self.overwrite, verbose=self.verbose)

        if not os.name == 'nt':
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
            self.framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        self.framework.defineElegantCommand(['elegant'])
        # if os.name == 'nt':
        #     self.framework.defineElegantCommand(['mpiexec','-np','10','pelegant'])
        self.framework.loadSettings(self.lattice_file)
        self.framework.change_Lattice_Code('All','elegant')

        ### Define starting lattice
        if self.start_lattice is None:
            if 'POSTINJ' in self.framework.latticeObjects:
                self.start_lattice = 'POSTINJ'
            else:
                self.start_lattice = self.framework.latticeObjects[0].objectname

        ### Apply any pre-tracking changes to elements
        if self.changes is not None:
            if isinstance(self.changes, (tuple, list)):
                for c in self.changes:
                    self.framework.load_changes_file(c)
            else:
                self.framework.load_changes_file(self.changes)

        ### Apply input arguments to element definitions
        ''' Apply arguments: [[element, parameter, value], [...]] '''
        for e, p, v in self.input_parameters:
            self.framework.modifyElement(e, p, v)

        ### Save the changes to the run directory
        self.framework.save_changes_file(filename=self.framework.subdirectory+'/changes.yaml', elements=self.input_parameters)

    def calculateBeamParameters(self):
        ### Are we running post-injector?
        if self.post_injector:
            ### Have we defined where the base files are?
            if self.base_files is not None:
                # print('Using base_files = ', self.base_files)
                self.framework[self.start_lattice].prefix = self.base_files
            ### If not, use the defaults based on the location of the CLARA example directory
            else:
                # print('Using CLARA_dir base_files = ', self.CLARA_dir+'/../../basefiles_'+str(self.scaling)+'/')
                self.framework[self.start_lattice].prefix = self.CLARA_dir+'/basefiles_'+str(self.scaling)+'/'

            ### Are we are setting the charge?
            if self.startcharge is not None:
                self.framework[self.start_lattice].bunch_charge = 1e-12*self.startcharge

            ### Are we are sub-sampling the distribution?
            if self.sample_interval is not None:
                # print('Sampling at ', self.sample_interval)
                self.framework[self.start_lattice].sample_interval = self.sample_interval#2**(3*4)

            ### TRACKING
            self.framework.track(startfile=self.start_lattice)
        else:
            self.framework.track()#startfile='FMS')


def optfunc(inputargs, dir=None, *args, **kwargs):
    global bestfit
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(inputargs, tmpdir, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(inputargs, dir, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return fitvalue
