import os
import sys
sys.path.append('../../../')
from SimulationFramework.Framework import *
import string
import re
import glob

cavity_regex = re.compile(r'(\d)\s+([-+]?\d+.\d+)+\s+([-+]?\d+.\d+)+\s+')

def get_Cavity_Phase(logfile):
    cavity_phases = []
    cavity_energy_gain = []
    with open(logfile) as file:
        n = 0
        for line in file:
            ans = cavity_regex.findall(line)
            if len(ans) > 0:
                cavity_phases.append(float(ans[-1][-1]))
                cavity_energy_gain.append(float(ans[-1][1]))
                n+=1
    return cavity_phases, cavity_energy_gain

def get_Cavity_Phases(settings, dir='.'):
    rcplattice = Framework('.', clean=False, verbose=False, overwrite=False)
    rcplattice.loadSettings(settings)
    files = filter(os.path.isfile, glob.glob(dir+"/*.log"))
    files.sort(key=lambda x: os.path.getmtime(x))
    cavity_phases = []
    cavity_energy_gain = []
    for f in files:
        section = os.path.basename(f).replace('.log','')
        # print 'section = ', section
        cavs = rcplattice[section].getElementType('cavity', setting='objectName')
        # phases = rcplattice['injector400'].getElementType('cavity', setting='phase')
        with open(f) as file:
            n = 0
            for line in file:
                ans = cavity_regex.findall(line)
                if len(ans) > 0:
                    cavity_phases.append([cavs[n], float(ans[-1][-1])])
                    cavity_energy_gain.append([cavs[n], float(ans[-1][1])])
                    n+=1
    return cavity_phases, cavity_energy_gain
