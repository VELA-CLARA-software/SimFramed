from ASTRAInjector import *
import blackbox as bb
# import sdds
import numpy as np
from constraints import *
import tempfile
import os
import shutil
import read_astra_file as raf

class fitnessFunc():

    def __init__(self, args):
        linac1field, linac1phase, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase = args
        print 'args = ', args
        self.parameters = args
        self.linacfields = [linac1field, linac2field, linac3field, linac4field]
        self.tmpdir = tempfile.mkdtemp(dir=os.getcwd())
        self.dirname = os.path.basename(self.tmpdir)
        astra = ASTRAInjector(self.dirname, overwrite=True)
        if not os.name == 'nt':
            astra.defineASTRACommand(['mpiexec','-np','12','/opt/ASTRA/astra_MPICH2.sh'])
        astra.loadSettings('short_240.settings')
        astra.modifySetting('linac1_field', linac1field)
        astra.modifySetting('linac1_phase', linac1phase)
        astra.modifySetting('linac2_field', linac2field)
        astra.modifySetting('linac2_phase', linac2phase)
        astra.modifySetting('linac3_field', linac3field)
        astra.modifySetting('linac3_phase', linac3phase)
        astra.modifySetting('4hc_field', fhcfield)
        astra.modifySetting('4hc_phase', fhcphase)
        astra.modifySetting('linac4_field', linac4field)
        astra.modifySetting('linac4_phase', linac4phase)
        astra.createInitialDistribution(npart=100, charge=250)
        astra.applySettings()
        astra.runASTRAFiles()
        ft = feltools(self.dirname)
        sddsfile = ft.convertToSDDS('test.in.128.4929.128')

    def removeTempDirectory(self):
        shutil.rmtree(self.tmpdir)

    def loadBeamFile(self):
        sddsbeam = sdds.SDDS(0)
        sddsbeam.load(self.dirname+'/test.in.128.4929.128.sdds')
        beam = {}
        for col in range(len(sddsbeam.columnName)):
            if len(sddsbeam.columnData[col]) == 1:
                beam[sddsbeam.columnName[col]] = sddsbeam.columnData[col][0]
            else:
                beam[sddsbeam.columnName[col]] = sddsbeam.columnData[col]
            SDDSparameterNames = list()
            # print beam
            for param in sddsbeam.columnName:
                if param in beam and isinstance(beam[param][0], (float, long)):
                    SDDSparameterNames.append(param)
        return beam

    def calculateBeamParameters(self):
        beam = raf.read_astra_beam_file(self.dirname+'/test.in.128.4929.128')
        c = [0] * 5
        w = [2,1,2,5,5]
        beam = self.loadBeamFile()
        sigmat = 1e12*np.std(beam['t'])
        sigmap = np.std(beam['p'])
        meanp = np.mean(beam['p'])
        fitp = 100*sigmap/meanp
        fhcfield = self.parameters[6]
        c[0] = 0 if sigmat < 0.5 else (np.abs(sigmat-0.5))
        c[1] = 0 if fitp < 0.5 else (np.abs(fitp-0.5))
        c[2] = 0 if meanp > 200 else (np.abs(meanp-200))
        c[3] = 0 if max(self.linacfields) < 32 else (np.abs(max(self.linacfields)-32))
        c[4] = 0 if fhcfield < 35 else (np.abs(fhcfield-35))
        fitness = map(lambda x,y: (x*y)**2,c,w)
        return np.sqrt(np.sum(c))

def optfunc(args):
    fit = fitnessFunc(args)
    fitvalue = fit.calculateBeamParameters()
    fit.removeTempDirectory()
    return fitvalue

# optfunc([21,20,25,-25,25,-25,30,187,25,0])

def main():
    bb.search(f=optfunc,  # given function
              box=[[10, 32], [-30,30], [10, 32], [-30,30], [10, 32], [-30,30], [10, 32], [135,200], [10, 32], [-30,30]],  # range of values for each parameter
              n=200,  # number of function calls on initial stage (global search)
              m=200,  # number of function calls on subsequent stage (local search)
              batch=4,  # number of calls that will be evaluated in parallel
              resfile='output.csv')  # text file where results will be saved

if __name__ == '__main__':
    main()
