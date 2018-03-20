import os
import numpy as np
from epics import caput
import math as m
import h5py

def varianceMLE(data, mean, value):
    s=0
    N=0
    for i in range(1, data[1].size):
        if data[9][i] > 0:
            # gets rid of bad particles
            s = s + m.pow(data[value][i] - mean, 2)
            N = N + 1
    return s / (N - 1)


def setBPM(pvRoot, data):
    print('    Setting BPM: ' + pvRoot)
    xArray = data[0]
    yArray = data[1]
    xAverage = np.average(xArray)
    yAverage = np.average(yArray)
    caput(pvRoot + ':X', xAverage)
    caput(pvRoot + ':Y', yAverage)


def setCamera(pvRoot, data):
    print('    Setting Camera: ' + pvRoot)
    xArray = data[0]
    yArray = data[1]
    xAverage = np.average(xArray)
    yAverage = np.average(yArray)
    xVar = varianceMLE(data, xAverage, 0)
    yVar = varianceMLE(data, yAverage, 1)
    caput(pvRoot + ':X', xAverage)
    caput(pvRoot + ':Y', yAverage)
    caput(pvRoot + ':SigmaX', m.sqrt(xVar))
    caput(pvRoot + ':SigmaY', m.sqrt(yVar))
    caput(pvRoot + ':DistribX', data[0])
    caput(pvRoot + ':DistribY', data[1])


def setWCM(pv, data):
    print('    Setting WCM: '+ pv)
    charge = 0
    for i in range(data[1].size):
        if data[9][i]>0:# gets rid of bad particles
            charge = charge + data[7][i]
    caput(pv, charge)


def setOutput(framework):
    for fileName in os.listdir(framework.subdir):
        if fileName.split('.')[-1] == 'hdf5':
            # load data
            f = h5py.File(fullName, 'r')
                listData = list(f['data'])
                f.close()
                data = np.array(image)
            data = np.loadtxt(fileName, unpack=True)
            pvRoot = framework.elements[fileName]['PV']
            if 'BPM' in fileName:
                setBPM(pvRoot, data)
            elif 'YAG' in fileName or 'SCR' in fileName:
                pvRoot = framework.elements[fileName]['Camera_PV']
                setCamera(pvRoot, data)
            elif 'SCOPE' in fileName:
                pv = pvRoot
                setWCM(pv, data)
            else:
                print 'Unregcognised output file.'
