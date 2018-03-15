import csv
from ASTRAInjector import *
import argparse

parser = argparse.ArgumentParser(description='Analyse ASTRA Optimisation Output')
parser.add_argument('-d', '--directory', default='twiss_best')
parser.add_argument('-s', '--settings', default='short_240_12b3.settings')


import read_beam_file as raf

beam = raf.beam()

def main():
    global app
    args = parser.parse_args()
    results = []

    with open(args.directory+'/twiss_best_solutions.csv', 'r') as csvfile:
      reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
      for row in reader:
        results.append(row)

    astra = ASTRAInjector('', overwrite=False)
    astra.loadSettings('short_240_12b3.settings')
    astra.fileSettings['test.2']['quad_K'] = results[0][0:6]
    astra.fileSettings['test.3']['quad_K'] = results[0][6:14]
    astra.fileSettings['test.4']['quad_K'] = results[0][14:16]
    astra.fileSettings['test.5']['quad_K'] = results[0][16:]

    astra.saveSettings(filename='short_240_12b3.settings')

    exit()

if __name__ == '__main__':
   main()




def runCSRTrackFiles( filename):
    command = csrTrackCommand + [filename]
    with open(os.devnull, "w") as f:
        subprocess.call(command, stdout=f, cwd=dirname)


dirname = 'twiss_best'

astra = ASTRAInjector(dirname, overwrite=True)
astra.loadSettings('short_240_12b3.settings')

ncpu = 40

if not os.name == 'nt':
    astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
    csrTrackCommand = ['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh']

astra.createInitialDistribution(npart=50000, charge=250)
if not os.name == 'nt':
    astra.fileSettings['test.5']['starting_distribution'] = 'end.fmt2.astra'
astra.applySettings()
if not os.name == 'nt':
    astra.runASTRAFiles(files=['test.1','test.2','test.3','test.4'])
    try:
        copyfile('csrtrk.in', dirname+'/'+'csrtrk.in')
    except:
        pass
    runCSRTrackFiles('csrtrk.in')
    beam.convert_csrtrackfile_to_astrafile(dirname+'/'+'end.fmt2', dirname+'/'+'end.fmt2.astra')

    astra.runASTRAFiles(files=['test.5'])
else:
    astra.runASTRAFiles()
