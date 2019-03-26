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
        try:
            shutil.rmtree(self.name)
        except:
            pass

class fitnessFunc():

    def __init__(self, args, tempdir, id=None, scaling=5, overwrite=True, verbose=False, summary=False, post_injector=True, startcharge=None, basefiles=None):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.scaling = scaling
        self.tmpdir = tempdir
        self.verbose = verbose
        self.summary = summary
        self.overwrite = overwrite
        self.post_injector = post_injector
        self.id = id
        self.startcharge = startcharge
        self.basefiles = basefiles
        print 'basefiles defined as ', self.basefiles
        ''' if only post-injector optimisation'''
        if self.post_injector:
            if len(args) == 10:
                dcp_factor = args[9]
            else:
                dcp_factor = 1
            linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase, bcangle = args[:9]
            self.parameters = dict(zip(['linac2field', 'linac2phase', 'linac3field', 'linac3phase', 'fhcfield', 'fhcphase', 'linac4field', 'linac4phase', 'bcangle'], args))
        else:
            ''' including injector parameters '''
            if len(args) == 16:
                dcp_factor = args[15]
            else:
                dcp_factor = 1
            gunphase, gunsol, linac1field, linac1phase, linac1sol1, linac1sol2, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase, bcangle = args[:15]
            self.parameters = dict(zip(['gunphase','gunsol','linac1field','linac1phase', 'linac1sol1', 'linac1sol2', 'linac2field', 'linac2phase', 'linac3field', 'linac3phase', 'fhcfield', 'fhcphase', 'linac4field', 'linac4phase', 'bcangle'], args))
        self.parameters['dcp_factor'] = dcp_factor
        self.npart=2**(3*scaling)
        ncpu = scaling*3
        if self.post_injector:
            self.sbandlinacfields = np.array([linac2field, linac3field, linac4field])
        else:
            self.sbandlinacfields = np.array([linac1field, linac2field, linac3field, linac4field])
        self.dirname = self.tmpdir
        self.framework = fw.Framework(self.dirname, overwrite=overwrite, verbose=verbose)
        if not os.name == 'nt':
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator'])
            self.framework.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        self.framework.defineElegantCommand(['elegant'])
        self.framework.loadSettings('Lattices/clara400_v12_v3_elegant_jkj.def')
        self.framework.load_changes_file('../Elegant/best_changes.yaml')
        if not self.post_injector:
            self.framework.generator.particles = self.npart
            self.framework.modifyElement('CLA-HRG1-GUN-CAV', 'phase', gunphase)
            self.framework.modifyElement('CLA-HRG1-GUN-SOL', 'field_amplitude', gunsol)
            self.framework.modifyElement('CLA-L01-CAV', 'field_amplitude', abs(linac1field))
            self.framework.modifyElement('CLA-L01-CAV', 'phase', linac1phase)
            self.framework.modifyElement('CLA-L01-CAV-SOL-01', 'field_amplitude', linac1sol1)
            self.framework.modifyElement('CLA-L01-CAV-SOL-02', 'field_amplitude', linac1sol2)
        self.framework.modifyElement('CLA-L02-CAV', 'field_amplitude', abs(linac2field))
        self.framework.modifyElement('CLA-L02-CAV', 'phase', linac2phase)
        self.framework.modifyElement('CLA-L03-CAV', 'field_amplitude', abs(linac3field))
        self.framework.modifyElement('CLA-L03-CAV', 'phase', linac3phase)
        self.framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', abs(fhcfield))
        self.framework.modifyElement('CLA-L4H-CAV', 'phase', fhcphase)
        self.framework.modifyElement('CLA-L04-CAV', 'field_amplitude', abs(linac4field))
        self.framework.modifyElement('CLA-L04-CAV', 'phase', linac4phase)
        self.framework['bunch_compressor'].set_angle(bcangle)
        self.framework.modifyElement('CLA-S07-DCP-01', 'factor', dcp_factor)
        self.framework.save_changes_file(filename=self.framework.subdirectory+'/changes.yaml')

    def calculateBeamParameters(self):
        # try:
            if self.post_injector:
                if self.basefiles is not None:
                    print 'basefiles defined as ', self.basefiles
                    self.framework['POSTINJ'].prefix = self.basefiles
                else:
                    print 'basefiles not defined! Using default location'
                    self.framework['POSTINJ'].prefix = '../../../basefiles_'+str(self.scaling)+'/'
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
