from ASTRAInjector import *
from CSRTrack import *
import numpy as np
from constraints import *
import os
import read_twiss_file as rtf
import read_beam_file as rbf
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import operator
import random
import multiprocessing
from scoop import futures
import csv

import shutil
import uuid

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

    def __init__(self, args, tempdir, npart=10000, ncpu=6, overwrite=True, verbose=False, summary=False):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.tmpdir = tempdir
        self.verbose = verbose
        self.summary = summary
        linac1field, linac1phase, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase, bcangle = args
        self.parameters = dict(zip(['linac1field','linac1phase', 'linac2field', 'linac2phase', 'linac3field', 'linac3phase', 'fhcfield', 'fhcphase', 'linac4field', 'linac4phase', 'bcangle'], args))
        self.npart = npart
        self.linacfields = [linac1field, linac2field, linac3field, linac4field]
        self.dirname = os.path.basename(self.tmpdir)
        self.astra = ASTRAInjector(self.dirname, overwrite=overwrite)
        self.csrtrack = CSRTrack(self.dirname, overwrite=overwrite)
        if not os.name == 'nt':
            self.astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.csrtrack.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        else:
            self.astra.defineASTRACommand(['astra'])
            self.csrtrack.defineCSRTrackCommand(['CSRtrack_1.201.wic.exe'])
        self.astra.loadSettings('short_240_12b3.settings')
        self.astra.modifySetting('linac4_field', abs(linac4field))
        self.astra.modifySetting('linac4_phase', linac4phase)
        self.astra.fileSettings['vb']['starting_distribution'] = '../CSRTest_Base/test.4.2573.001'
        self.astra.fileSettings['vb']['variable_bunch_compressor']['angle'] = abs(bcangle)
        self.astra.fileSettings['vb']['LSPCH'] = True
        self.astra.fileSettings['vb']['LSPCH3D'] = True

    def between(self, value, minvalue, maxvalue, absolute=True):
        if absolute:
            result = max([minvalue,min([maxvalue,abs(value)])])
        else:
            result = np.sign(value)*max([minvalue,min([maxvalue,abs(value)])])
        return result

    def calculateBeamParameters(self):
        bcangle = float(self.astra.fileSettings['vb']['variable_bunch_compressor']['angle'])
        try:
            if bcangle < 0.05 or bcangle > 0.12:
                raise ValueError
            self.astra.createInitialDistribution(npart=self.npart, charge=250)
            ''' Modify the last file to use to CSRTrack output as input'''
            self.astra.fileSettings['test.5']['starting_distribution'] = 'end.fmt2.astra'
            self.astra.applySettings()
            ''' Run ASTRA upto VBC '''
            self.astra.runASTRAFiles(['vb'])
            ''' Write Out the CSRTrack file based on the BC angle (assumed to be 0.105) '''
            self.csrtrack.writeCSRTrackFile('csrtrk.in', angle=bcangle, forces='projected', inputfile='../CSRTest_Base/test.4.2573.001')
            ''' Run CSRTrack'''
            self.csrtrack.runCSRTrackFile('csrtrk.in')
            ''' Convert CSRTrack output file back in to ASTRA format '''
            self.beam.convert_csrtrackfile_to_astrafile(self.dirname+'/'+'end.fmt2', self.dirname+'/'+'end.fmt2.astra')
            ''' Run the next section of the lattice in ASTRA, using the CSRTrack output as input '''
            self.astra.runASTRAFiles(files=['test.5'])

            self.beam.read_astra_beam_file(self.dirname+'/test.5.4936.001')
            self.beam.slice_length = 0.03e-12
            self.beam.bin_time()
            sigmat = 1e12*np.std(self.beam.t)
            sigmap = np.std(self.beam.p)
            meanp = np.mean(self.beam.p)
            emitx = 1e6*self.beam.normalized_horizontal_emittance
            emity = 1e6*self.beam.normalized_vertical_emittance
            fitp = 100*sigmap/meanp
            fhcfield = self.parameters['fhcfield']
            peakI, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum = self.beam.sliceAnalysis()
            constraintsList = {
                'peakI': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 400, 'weight': 100},
                'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.05, 'weight': 3},
                'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 0.75, 'weight': 2.5},
                'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 0.75, 'weight': 2.5},
                'peakIMomentum': {'type': 'greaterthan','value': 1e-6*peakIMomentum, 'limit': 240, 'weight': 10},
                'linac fields': {'type': 'lessthan', 'value': self.linacfields, 'limit': 32, 'weight': 100},
                '4hc field': {'type': 'lessthan', 'value': fhcfield, 'limit': 35, 'weight': 100},
                'horizontal emittance': {'type': 'lessthan', 'value': emitx, 'limit': 1, 'weight': 4},
                'vertical emittance': {'type': 'lessthan', 'value': emity, 'limit': 1, 'weight': 4},
            }
            self.twiss.read_astra_emit_files(self.dirname+'/test.5.Zemit.001')
            constraintsList5 = {
                'last_exn_5': {'type': 'lessthan', 'value': 1e6*self.twiss['enx'], 'limit': 1, 'weight': 1},
                'last_eyn_5': {'type': 'lessthan', 'value': 1e6*self.twiss['eny'], 'limit': 1, 'weight': 1},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsList5)
            fitness = self.cons.constraints(constraintsList)
            if self.verbose:
                print self.cons.constraintsList(constraintsList)
            if self.summary:
                self.astra.createHDF5Summary(reference='CSR')
            return fitness
        except:
            return 1e6

def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return (fitvalue,)

# results = []
#
# with open('longitudinal_best/longitudinal_best_solutions.csv', 'r') as csvfile:
#     reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
#     for row in reader:
#         results.append(row)
# best = results[0]
astra = ASTRAInjector('twiss_best', overwrite=False)
astra.loadSettings('short_240_12b3.settings')
parameters = []
parameters.append(astra.getSetting('linac1_field')[0][1])
parameters.append(astra.getSetting('linac1_phase')[0][1])
parameters.append(astra.getSetting('linac2_field')[0][1])
parameters.append(astra.getSetting('linac2_phase')[0][1])
parameters.append(astra.getSetting('linac3_field')[0][1])
parameters.append(astra.getSetting('linac3_phase')[0][1])
parameters.append(astra.getSetting('4hc_field')[0][1])
parameters.append(astra.getSetting('4hc_phase')[0][1])
parameters.append(astra.getSetting('linac4_field')[0][1])
parameters.append(astra.getSetting('linac4_phase')[0][1])
parameters.append(astra.fileSettings['vb']['variable_bunch_compressor']['angle'])
best = parameters
print optfunc(best, dir=os.getcwd()+'/longitudinal_best_CSR', npart=2**16, ncpu=20, overwrite=True, verbose=True, summary=True)
exit()
