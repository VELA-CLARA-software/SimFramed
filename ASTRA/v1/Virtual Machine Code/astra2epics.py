import sys,os
import numpy as np
os.system('PATH=$PATH:/home/vsim/Desktop/epics/base/bin/linux-x86_64')
from epics import caget,caput
os.environ["EPICS_CA_SERVER_PORT"]="6000"
import math as m

def varianceMLE(data,mean,value):
	s=0
	N=0
	for i in range(1,data[1].size):
		if data[9][i]>0:# gets rid of bad particles
			s = s + m.pow(data[value][i]-mean,2)
			N = N + 1
	return s/(N-1)

def setBPM(pvRoot,data):
	print('	Setting BPM: '+pvRoot)
	xArray = data[0]
	yArray = data[1]
	xAverage = np.average(xArray)
	yAverage = np.average(yArray)
	caput(pvRoot+':X',xAverage)
	caput(pvRoot+':Y',yAverage)

def setCamera(pvRoot,data):
	print('	Setting Camera: '+pvRoot)
	xArray = data[0]
	yArray = data[1]
	xAverage = np.average(xArray)
	yAverage = np.average(yArray)
	xVar = varianceMLE(data,xAverage,0)
	yVar = varianceMLE(data,yAverage,1)
	caput(pvRoot+':X',xAverage)
	caput(pvRoot+':Y',yAverage)
	caput(pvRoot+':SigmaX',m.sqrt(xVar))
	caput(pvRoot+':SigmaY',m.sqrt(yVar))
	caput(pvRoot+':DistribX',data[0])
	caput(pvRoot+':DistribY',data[1])

def setWCM(pv,data):
	print('	Setting WCM: '+pv)
	charge=0
	for i in range(data[1].size):
		if data[9][i]>0:# gets rid of bad particles
			charge = charge + data[7][i]	
	caput(pv,charge)


def setOutput(inFileName):
	dir_path = os.path.dirname(os.path.realpath(__file__))
	for fileName in os.listdir(dir_path):
		if fileName.startswith(inFileName[:-3]) and fileName.split('.')[1][:2]=='VM':
			#load data
			data = np.loadtxt(dir_path+'/'+fileName,unpack=True)
			pvRoot = fileName.split('.')[1]
			if 'BPM' in fileName:
				setBPM(pvRoot,data)
			elif 'CAM' in fileName:
				setCamera(pvRoot,data)
			elif 'SCOPE' in fileName:
				pv=pvRoot
				setWCM(pv,data)

	
