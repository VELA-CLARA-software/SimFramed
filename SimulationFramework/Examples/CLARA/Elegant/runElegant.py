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
        self._lattice_file = lattice
        self._base_files = None
        self._verbose = False

    def set_CLARA_directory(self, clara_dir):
        self.CLARA_dir = clara_dir

    @property
    def base_files(self):
        return self._base_files
    @base_files.setter
    def base_files(self, location):
        self._base_files = location

    @property
    def lattice_file(self):
        return self._lattice_file
    @lattice_file.setter
    def lattice_file(self, file):
        self._lattice_file = file

    @property
    def verbose(self):
        return self._verbose
    @verbose.setter
    def verbose(self, file):
        self._verbose = file

    def setup_lattice(self, inputargs, tempdir, scaling=5, overwrite=True, post_injector=True, startcharge=None, basefiles=None, changes=None, verbose=False, sample_interval=1, *args, **kwargs):
        if 'lattice' in kwargs:
            # print 'loading lattice ', kwargs['lattice']
            self.lattice_file = kwargs['lattice']
        self.scaling = scaling
        self.dirname = tempdir
        self.overwrite = overwrite
        self.post_injector = post_injector
        self.startcharge = startcharge
        self.sample_interval = sample_interval
        if basefiles is not None:
            self.base_files = basefiles
        self.input_parameters = inputargs

        self.npart=2**(3*scaling)
        ncpu = scaling*3

        self.framework = fw.Framework(self.dirname, overwrite=overwrite, verbose=verbose)

        if not os.name == 'nt':
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
            self.framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        self.framework.defineElegantCommand(['elegant'])

        self.framework.loadSettings(self.lattice_file)
        self.framework.change_Lattice_Code('All','elegant')
        if 'POSTINJ' in self.framework.latticeObjects:
            self.start_lattice = 'POSTINJ'
        else:
            self.start_lattice = 'S02'
        if changes is not None:
            if isinstance(changes, (tuple, list)):
                for c in changes:
                    # print 'loading changes file: ', c
                    self.framework.load_changes_file(c)
            else:
                self.framework.load_changes_file(changes)

        ''' Apply arguments: [[element, parameter, value], [...]] '''
        # self.input_parameter_names = []
        for e, p, v in self.input_parameters:
            # self.input_parameter_names.append(e)
            # if e == 'bunch_compressor' and p == 'set_angle':
            #     print 'modifying VBC = ', v
            #     self.framework['bunch_compressor'].set_angle(float(v))
            # else:
            # print 'modifying ',e,'[',p,'] = ', v
            self.framework.modifyElement(e, p, v)

        # pnames = []
        # if not hasattr(self, 'parameter_names'):
        #     for e in self.input_parameter_names:
        #         if e == 'bunch_compressor':
        #             pnames.append('CLA-VBC-MAG-DIP-01')
        #             pnames.append('CLA-VBC-MAG-DIP-02')
        #             pnames.append('CLA-VBC-MAG-DIP-03')
        #             pnames.append('CLA-VBC-MAG-DIP-04')
        #         else:
        #             pnames.append(e)
        # else:
        #     pnames= self.parameter_names
        self.framework.save_changes_file(filename=self.framework.subdirectory+'/changes.yaml', elements=self.input_parameters)

    def calculateBeamParameters(self):
        if self.post_injector:
            if self.base_files is not None:
                self.framework[self.start_lattice].prefix = self.base_files
            else:
                self.framework[self.start_lattice].prefix = self.CLARA_dir+'/../../basefiles_'+str(self.scaling)+'/'
            if self.startcharge is not None:
                self.framework[self.start_lattice].bunch_charge = 1e-12*self.startcharge
            if self.sample_interval is not None:
                self.framework[self.start_lattice].sample_interval = self.sample_interval#2**(3*4)

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
