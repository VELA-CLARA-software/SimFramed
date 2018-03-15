from ASTRAInjector import *
import argparse

parser = argparse.ArgumentParser(description='Create HDF5 Summary file from a directory')
parser.add_argument('directory')
parser.add_argument('-r', '--reference', default='')
parser.add_argument('-np', '--numpart', default=1000)
parser.add_argument('-q', '--charge', default=250)
parser.add_argument('-s', '--settings', default='short_240_12b3.settings')

def createHDF5File(directory, settings, reference, numpart, charge):

    astra = ASTRAInjector(directory, overwrite=False)
    astra.loadSettings(settings)
    astra.globalSettings['npart'] = int(numpart)
    astra.globalSettings['charge'] = int(charge)/1000.0
    astra.createHDF5Summary(reference=reference)

def main():
    args = parser.parse_args()
    createHDF5File(directory=args.directory, reference=args.reference, settings=args.settings, numpart=args.numpart, charge=args.charge)

if __name__ == '__main__':
   main()
