
import sys,os
import editOutputFiles
import astra2epics
import time

from epics import caget,caput


def find(string, infile):
	for i,line in enumerate(infile):
		if string in line:
			return i

def setZOffset(fileName):
	import numpy as np
	data = np.loadtxt(dir_path+'/temp-start.ini',unpack=True).T
	original = open(dir_path+'/temp-'+fileName, 'r')
	lines = original.read().split("\n")
	original.close()
	place = find('ZSTART=',lines)
	offset=lines[place].split('=')[1]
	print('Offset in next sim will be = '+str(float(offset)))
	data[0][2]=float(offset)
	np.savetxt(dir_path+'/temp-start.ini',data, delimiter='  ', newline='\n')


inList = sys.argv[1].split(',')

#Set old values back to zero        not elegent but works :)
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)
for fileName in os.listdir(dir_path):
	#find files used to set epics PVs
	if '.VM' in fileName:
		pvRoot = fileName.split('.')[1]
		#Reset BPMs
		if 'BPM' in fileName:
			print(fileName)
			caput(pvRoot+':X',0.0)
			caput(pvRoot+':Y',0.0)
		#Reset Camera
		if 'CAM' in fileName:
			print(fileName)
			#reset that BPM data
			caput(pvRoot+':X',0.0)
			caput(pvRoot+':Y',0.0)
			caput(pvRoot+':SigmaX',0.0)
			caput(pvRoot+':SigmaY',0.0)
			caput(pvRoot+':DistribX',[])
			caput(pvRoot+':DistribY',[])

for fileName in os.listdir(dir_path):
	if fileName.split('.')[-1].isdigit() or fileName.split('.')[0]=='C1' or fileName.split('.')[0]=='V2' or fileName.split('.')[0]=='CV' or fileName.split('.')[0]=='V1' or fileName.split('.')[0]=='SP':
		os.system('rm '+dir_path+'/'+fileName)

for index,section in enumerate(inList):
	#run ASTRA for a section
	print(section)
	print('Running ASTRA for section '+section[:-3]+'...')
	os.system(dir_path+'/Astra  /home/vmsim/Desktop/V2/ASTRA/temp-'+str(section)+' > /home/vmsim/Desktop/V2/ASTRA/temp-'+str(section[:-3])+'.out')
	original = open('/home/vmsim/Desktop/V2/ASTRA/temp-'+str(section[:-3])+'.out', 'r')
	lines = original.read().split("\n")
	original.close()
	for i in range(5,2,-1):
		print(lines[-i])
	print('Renaming files...')
	editOutputFiles.rename(str(section))
	print('Setting EPICS values with ASTRA output...')
	astra2epics.setOutput(str(section))

	

	dir_path = os.path.dirname(os.path.realpath(__file__))
	for fileName in os.listdir(dir_path):
		if fileName.startswith('temp-'+section[:-3]) and fileName.split('.')[1].isdigit():
			if section!=inList[-1] and section=='CV.in' and inList[index+1]=='SP.in':
				#going straight through therefore dont' need to rotate
				os.system('cp '+dir_path+'/'+fileName+' '+dir_path+'/temp-start.ini')
				setZOffset(inList[index+1])
			elif section!=inList[-1] and section=='C1.in' and inList[index+1]=='C2.in':
				#going straight through therefore dont' need to rotate
				os.system('cp '+dir_path+'/'+fileName+' '+dir_path+'/temp-start.ini')
				setZOffset(inList[index+1])
			elif section!=inList[-1] and section=='V1.in' and inList[index+1]=='V2.in':
				#going straight through therefore dont' need to rotate
				os.system('cp '+dir_path+'/'+fileName+' '+dir_path+'/temp-start.ini')
				setZOffset(inList[index+1])
			elif section!=inList[-1] and section=='V1.in' and inList[index+1]=='SP.in':
				editOutputFiles.rotateAndSave(fileName,'R')
				setZOffset(inList[index+1])
			elif section!=inList[-1] and section=='C1.in' and inList[index+1]=='CV.in':
				editOutputFiles.rotateAndSave(fileName,'R')
				setZOffset(inList[index+1])
			elif section!=inList[-1] and section=='CV.in' and inList[index+1]=='V2.in':
				editOutputFiles.rotateAndSave(fileName,'L')
				setZOffset(inList[index+1])
			else:
				print('Junction doesn\'t exsit or this is the last file')


os.system('cp /home/vmsim/Desktop/V2/ASTRA/temp-start_BACKUP.ini /home/vmsim/Desktop/V2/ASTRA/temp-start.ini')
