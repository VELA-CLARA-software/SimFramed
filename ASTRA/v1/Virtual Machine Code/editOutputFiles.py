import os,sys
import numpy as np
import math as m
import shutil


def rename(inFileName):
	dir_path = os.path.dirname(os.path.realpath(__file__))
	#print('Looking in: '+dir_path)
	fIn = open(dir_path+'/temp-'+inFileName, 'r')
	inLines = fIn.read().split("\n")
	fIn.close()
	#loop throught the output files
	for fileName in os.listdir(dir_path):
		#check if file is from the right section oof the line and it is a screen output file
   		if fileName.startswith('temp-'+inFileName[:-3]) and fileName.split('.')[1].isdigit():
			#print(fileName)
			data = np.loadtxt(dir_path+'/'+fileName,unpack=True).T
			zpos = round(data[0][2],2)
			#zpos=float(fileName.split('.')[1])/100
			#print('ah')
			#print(zpos)
			for i,line in enumerate(inLines):
				if 'Screen(' in line:
					#print(round(float(line.split('=')[1]),2))
					if round(float(line.split('=')[1]),2)==zpos:
						elementName = inLines[i-2].split('!')[1]
						os.system('mv '+dir_path+'/'+fileName+' '+dir_path+'/'+inFileName[:-2]+str(elementName[:-1]))
						print('	Moving: '+fileName+' > '+inFileName[:-2]+str(elementName[:-1]))
					else:
						continue#print('Round read positionion were not the same.')
    		else:
        		continue

def rotate(zIN,xIN, D):
	#always rotate 45 degrees
	#coord[0]=z/pz,coord[1]=x/px
	const= m.sqrt(2)/2
	if D=='R':
		z=const*(zIN-xIN)
		x=const*(xIN+zIN)
	elif D=='L':
		z=const*(zIN+xIN)
		x=const*(xIN-zIN)
	else:
		print('Error: Not a valid direction for rotating 45 degrees.')
	return [z,x]

def rotateAndSave(fileName,D):
	dir_path = os.path.dirname(os.path.realpath(__file__))

	data = np.loadtxt(dir_path+'/'+fileName,unpack=True).T

	refpz = data[0][5]
	print('Reference Z mom : '+str(refpz))

	for i in range(len(data)):
		#catch condition if reference particle

		if i==0:
			#rotate z,x
			ZXRot = rotate(0,data[i][0],D)
			if D=='R':
				ZXRot[1] = ZXRot[1]-0.1775
			elif D=='L':
				ZXRot[1] = ZXRot[1]+0.1775
			#ZXRot[0]=-0.38845
			#rotate pz,px
			momZXRot = rotate(data[i][5],data[i][3],D)
			#momZXRot[0] = momZXRot[0]+data[i][5]
		else:
			#rotate z,x
			ZXRot = rotate(data[i][2],data[i][0],D)
			if D=='R':
				ZXRot[1] = ZXRot[1]-0.1775
			elif D=='L':
				ZXRot[1] = ZXRot[1]+0.1775
			#rotate pz,px
			momZXRot = rotate(data[i][5]+refpz,data[i][3],D)
			momZXRot[0] = momZXRot[0]-data[0][5]

		#Replace the old values with teh new rotated ones
		data[i] = [ZXRot[1],data[i][1],ZXRot[0],momZXRot[1],data[i][4],momZXRot[0],data[i][6],data[i][7],data[i][8],data[i][9]]

	np.savetxt(dir_path+'/temp-start.ini',data, delimiter='  ', newline='\n')
	np.savetxt(dir_path+'/dipoleStop.ini',data, delimiter='  ', newline='\n')
