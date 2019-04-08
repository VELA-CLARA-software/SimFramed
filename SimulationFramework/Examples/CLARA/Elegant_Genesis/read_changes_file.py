import numpy as np
import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
import SimulationFramework.Framework as fw
import argparse

parser = argparse.ArgumentParser(description='Print out a changes file')
parser.add_argument('filename', metavar='filename',
                   help='Changes file to read in', type=str)
parser.add_argument('-l', '--lattice', metavar='lattice',
                   help='Framework Lattice to load', type=str, default='Lattices/clara400_v12_v3.def')
parser.add_argument('-o', '--output', metavar='output',
                   help='Changes file to write out', type=str, default='')

args = parser.parse_args()

framework = fw.Framework(None, overwrite=False, verbose=False)
framework.loadSettings(args.lattice)
d = framework.load_changes_file(filename=args.filename, apply=True)
if args.output is not '':
        framework.save_changes_file(filename=args.output)
else:
    framework.save_changes_file(filename=args.filename+'.tmp')
