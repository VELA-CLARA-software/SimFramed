import os
import numpy as np
from epics import caput
import math as m

from SimulationFramework.Modules import read_beam_file as rbf


def varianceMLE(array, mean):
    s = 0
    N = 0
    for i in array:
        # gets rid of bad particles
        s = s + m.pow(i - mean, 2)
        N = N + 1
    return s / (N - 1)


def setBPM(pvRoot, beam):
    print('    Setting BPM: ' + pvRoot)
    xArray = beam.x
    yArray = beam.y
    xAverage = np.average(xArray)
    yAverage = np.average(yArray)
    caput(pvRoot + ':X', xAverage)
    caput(pvRoot + ':Y', yAverage)


def setCamera(pvRoot, beam):
    print('    Setting Camera: ' + pvRoot)
    xArray = beam.x
    yArray = beam.y
    xAverage = np.average(xArray)
    yAverage = np.average(yArray)
    xVar = varianceMLE(xArray, xAverage)
    yVar = varianceMLE(yArray, yAverage)
    caput(pvRoot + ':ANA:X_RBV', xAverage)
    caput(pvRoot + ':ANA:Y_RBV', yAverage)
    caput(pvRoot + ':ANA:SigmaX_RBV', m.sqrt(xVar))
    caput(pvRoot + ':ANA:SigmaY_RBV', m.sqrt(yVar))
    #caput(pvRoot + ':DistribX', xArray)
    #caput(pvRoot + ':DistribY', yArray)


def setWCM(pv, data):
    print('    Setting WCM: ' + pv)
    charge = 0
    for i in range(data[1].size):
        if data[9][i]>0:# gets rid of bad particles
            charge = charge + data[7][i]
    caput(pv, charge)


def setOutputs(framework):
    beam = rbf.beam()
    for fileName in os.listdir(framework.subdir):
        if fileName.split('.')[-1] == 'hdf5':
            # load data
            name = fileName.split('.')[0]
            #print name
            beam.read_HDF5_beam_file(framework.subdir + '/' + fileName,
                                     local=True)
            if 'BPM' in fileName:
                pvRoot = 'VM-' + framework.getElement(element=name,
                                                      setting='PV',
                                                      default='Null')
                setBPM(pvRoot, beam)
            elif 'YAG' in fileName or 'SCR' in fileName:
                #print fileName
                pvRoot = 'VM-' + framework.getElement(element=name,
                                                      setting='camera_PV',
                                                      default='Null')
                setCamera(pvRoot, beam)
            elif 'SCOPE' in fileName:
                pvRoot = 'VM-' + framework.getElement(element=name,
                                                      setting='PV',
                                                      default='Null')
                setWCM(pvRoot, beam)
            else:
                print '    Unregcognised output file', fileName
