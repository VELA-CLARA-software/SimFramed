import numpy as np
import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *
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

    def __init__(self, args, tempdir, scaling=4, overwrite=True, verbose=False, summary=False, post_injector=True):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.scaling = scaling
        self.tmpdir = tempdir
        self.overwrite = overwrite
        linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase, bcangle = args
        self.parameters = dict(zip(['linac2field', 'linac2phase', 'linac3field', 'linac3phase', 'fhcfield', 'fhcphase', 'linac4field', 'linac4phase', 'bcangle'], args))
        self.npart=2**(3*scaling)
        ncpu = scaling*3
        self.sbandlinacfields = np.array([linac2field, linac3field, linac4field])
        self.dirname = os.path.basename(self.tmpdir)
        self.framework = Framework(self.dirname, overwrite=overwrite, verbose=verbose)
        if not os.name == 'nt':
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
            self.framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        self.framework.defineElegantCommand(['elegant'])

        self.framework.loadSettings('Lattices/clara400_v12_v3_elegant.def')
        self.framework.modifyElement('CLA-L02-CAV', 'field_amplitude', abs(linac2field))
        self.framework.modifyElement('CLA-L02-CAV', 'phase', linac2phase)
        self.framework.modifyElement('CLA-L03-CAV', 'field_amplitude', abs(linac3field))
        self.framework.modifyElement('CLA-L03-CAV', 'phase', linac3phase)
        self.framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', abs(fhcfield))
        self.framework.modifyElement('CLA-L4H-CAV', 'phase', fhcphase)
        self.framework.modifyElement('CLA-L04-CAV', 'field_amplitude', abs(linac4field))
        self.framework.modifyElement('CLA-L04-CAV', 'phase', linac4phase)
        self.framework['bunch_compressor'].set_angle(abs(bcangle))


    def calculateBeamParameters(self):
        try:
            self.framework['POSTINJ'].file_block['input']['prefix'] = '../basefiles_'+str(self.scaling)+'/'
            self.framework.track()#startfile='FMS')
            self.beam.read_HDF5_beam_file(self.dirname+'/CLA-FMS-APER-01.hdf5')
            ## CONVERT THIS TO A DIST FILE!!!
        except Exception as e:
            print(e)
            return 1e6

def optfunc(inputargs, dir=None, *args, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(inputargs, tmpdir, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(inputargs, dir, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return fitvalue
