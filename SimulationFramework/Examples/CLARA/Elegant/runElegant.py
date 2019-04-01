import numpy as np
import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
import SimulationFramework.Framework as fw
from SimulationFramework.Modules.constraints import *
import SimulationFramework.Modules.read_beam_file as rbf
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

class fitnessFunc():

    def __init__(self):
        self.CLARA_dir = os.path.relpath(os.path.dirname(os.path.dirname( os.path.abspath(__file__))))# if CLARA_dir is None else CLARA_dir
        self.cons = constraintsClass()
        self.beam = rbf.beam()

    def set_CLARA_directory(self, clara_dir):
        self.CLARA_dir = clara_dir

    def setup_lattice(self, inputargs, tempdir, scaling=5, overwrite=True, verbose=False, post_injector=True, startcharge=None, basefiles=None, changes=None, *args, **kwargs):
        self.scaling = scaling
        self.dirname = tempdir
        self.verbose = verbose
        self.overwrite = overwrite
        self.post_injector = post_injector
        self.startcharge = startcharge
        self.basefiles = basefiles
        self.parameters = inputargs

        self.npart=2**(3*scaling)
        ncpu = scaling*3

        # if self.post_injector:
        #     self.sbandlinacfields = np.array([self.parameters[''], linac3field, linac4field])
        # else:
        #     self.sbandlinacfields = np.array([linac1field, linac2field, linac3field, linac4field])

        self.framework = fw.Framework(self.dirname, overwrite=overwrite, verbose=verbose)

        if not os.name == 'nt':
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
            self.framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        self.framework.defineElegantCommand(['elegant'])

        self.framework.loadSettings('Lattices/clara400_v12_v3_elegant_jkj.def')
        if changes is not None:
            self.framework.load_changes_file(changes)#self.CLARA_dir+'/Elegant/best_changes.yaml')

        ''' Apply arguments: [[element, parameter, value], [...]] '''
        for e, p, v in self.parameters:
            print e,p,v
            if e == 'bunch_compressor' and p == 'set_angle':
                self.framework['bunch_compressor'].set_angle(float(v))
            else:
                self.framework.modifyElement(e, p, v)

        self.framework.save_changes_file(filename=self.framework.subdirectory+'/changes.yaml', function=float)

    def calculateBeamParameters(self):
        # try:
            if self.post_injector:
                if self.basefiles is not None:
                    print 'basefiles defined as ', self.basefiles
                    self.framework['POSTINJ'].prefix = self.basefiles
                else:
                    print 'basefiles not defined! Using default location'

                    self.framework['POSTINJ'].prefix = self.CLARA_dir+'/../../basefiles_'+str(self.scaling)+'/'
                if self.startcharge is not None:
                    self.framework['POSTINJ'].bunch_charge = 1e-12*self.startcharge
                self.framework.track(startfile='POSTINJ')
            else:
                self.framework.track()#startfile='FMS')
            # self.beam.read_HDF5_beam_file(self.dirname+'/CLA-FMS-APER-01.hdf5')
        # except Exception as e:
        #     print(e)
        #     return 1e6

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
