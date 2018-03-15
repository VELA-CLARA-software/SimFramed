from ASTRAInjector import *
from CSRTrack import *
import numpy as np
from constraints import *
import os
import read_twiss_file as rtf
import read_beam_file as rbf
import csv

class fitnessFunc():

    def __init__(self,  tempdir, ncpu=6, overwrite=True, verbose=False, summary=False, csrtrack=True, csrforce='projected', vbsc=True, charge=250):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.tmpdir = tempdir
        self.verbose = verbose
        self.summary = summary
        self.npart = 2**16
        self.dirname = self.tmpdir
        self.astra = ASTRAInjector(self.dirname, overwrite=overwrite)
        self.csrtrack = CSRTrack(self.dirname, overwrite=overwrite)
        if not os.name == 'nt':
            self.astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            self.csrtrack.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        else:
            self.astra.defineASTRACommand(['astra'])
            self.csrtrack.defineCSRTrackCommand(['CSRtrack_1.201.wic.exe'])
        self.astra.loadSettings('short_240_12b3.settings')
        bcangle = self.astra.fileSettings['vb']['variable_bunch_compressor']['angle']
        self.astra.fileSettings['vb']['starting_distribution'] = '../CSRTest_Base/test.4.2573.001'
        self.astra.fileSettings['vb']['LSPCH'] = vbsc
        self.astra.fileSettings['vb']['LSPCH3D'] = vbsc
        try:
            self.astra.createInitialDistribution(npart=self.npart, charge=charge)
            if csrtrack:
                ''' Modify the last file to use to CSRTrack output as input'''
                self.astra.fileSettings['test.5']['starting_distribution'] = 'end.fmt2.astra'
            self.astra.applySettings()
            ''' Run ASTRA upto VBC '''
            self.astra.runASTRAFiles(['vb'])
            if csrtrack:
                ''' Write Out the CSRTrack file based on the BC angle (assumed to be 0.105) '''
                self.csrtrack.writeCSRTrackFile('csrtrk.in', angle=bcangle, forces=csrforce, inputfile='../CSRTest_Base/test.4.2573.001')
                ''' Run CSRTrack'''
                self.csrtrack.runCSRTrackFile('csrtrk.in')
                ''' Convert CSRTrack output file back in to ASTRA format '''
                self.beam.convert_csrtrackfile_to_astrafile(self.dirname+'/'+'end.fmt2', self.dirname+'/'+'end.fmt2.astra')
            ''' Run the next section of the lattice in ASTRA, using the CSRTrack output as input '''
            self.astra.runASTRAFiles(files=['test.5'])

            if self.summary:
                reference = tempdir
                if os.name == 'nt':
                    reference = reference + '_windows'
                self.astra.createHDF5Summary(reference=reference)
        except:
            print 'Error!'

''' ASTRA No CSR, No SC '''
# fitnessFunc('ASTRA_NoCSR_NoSC', ncpu=20, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=False)

''' CSRTrack No CSR, No SC '''
# fitnessFunc('CSRTrack_NoSC', ncpu=20, overwrite=True, verbose=False, summary=True, csrtrack=True, csrforce='none', vbsc=True)

''' ASTRA No CSR, 3D-SC '''
# fitnessFunc('ASTRA_NoCSR', ncpu=20, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=True)

''' ASTRA No CSR, 3D-SC Low Q'''
# fitnessFunc('ASTRA_NoCSR_10pc', ncpu=40, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=True, charge=10)
# fitnessFunc('ASTRA_NoCSR_1pc', ncpu=40, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=True, charge=1)
# fitnessFunc('ASTRA_NoCSR_100pc', ncpu=40, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=True, charge=100)
# fitnessFunc('ASTRA_NoCSR_250pc', ncpu=40, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=True, charge=250)
# fitnessFunc('ASTRA_NoCSR_250pc_20cores', ncpu=20, overwrite=True, verbose=False, summary=True, csrtrack=False, vbsc=True, charge=250)

''' CSRTrack CSR, No SC '''
# fitnessFunc('CSRTrack_CSR', ncpu=40, overwrite=True, verbose=False, summary=True, csrtrack=True, csrforce='projected', vbsc=True)

''' CSRTrack 3D-CSR, 3D-SC '''
# fitnessFunc('CSRTrack_CSR_SC', ncpu=20, overwrite=True, verbose=False, summary=True, csrtrack=True, csrforce='csr_g_to_p', vbsc=True)
