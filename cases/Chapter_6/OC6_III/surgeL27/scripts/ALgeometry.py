import re
import pandas as pd
from natsort import natsorted
import glob
import os
import sys

# FUNCTION DEFINITION
def getPosition(filename, index):
	"""
	This function opens filename,
    and returns position values at 
    position 'index'
    """
	df = pd.read_csv(filename).to_dict()
	x = df.get('x')
	y = df.get('y')
	z = df.get('z')
	x_value = list(x.values())[index]
	y_value = list(y.values())[index]
	z_value = list(z.values())[index]

	return x_value, y_value, z_value;


def transformPoints(x,y,z):
	"""
	Given the AL center coordinates,
	return the AL sides coordinates. 
    """

	n = len(x)
	if n>1:
		x_new = [0]*(n+1)
		y_new = [0]*(n+1)
		z_new = [0]*(n+1)

		x_new[0] = 1.5*x[0] - 0.5*x[1]
		y_new[0] = 1.5*y[0] - 0.5*y[1]
		z_new[0] = 1.5*z[0] - 0.5*z[1]

		x_new[-1] = 1.5*x[-1] - 0.5*x[-2]
		y_new[-1] = 1.5*y[-1] - 0.5*y[-2]
		z_new[-1] = 1.5*z[-1] - 0.5*z[-2]

		for i in range(1,n):
			x_new[i] = 0.5*x[i-1] + 0.5*x[i]
			y_new[i] = 0.5*y[i-1] + 0.5*y[i]
			z_new[i] = 0.5*z[i-1] + 0.5*z[i]	 	

	else:
		x_new = x
		y_new = y
		z_new = z

	return x_new, y_new, z_new

def getClosestTime(filename, tref):
	"""
	This function opens filename,
    and looks for the time value 
    closest to 'tref' 
    """
	df = pd.read_csv(filename).to_dict()
	t = df.get('time')
	# Get time value closest to ref_t
	index, time_value = min(t.items(), key=lambda x: abs(tref - x[1]))
	return index, time_value;


def getFilenames(regionName, componentName):
	"""
	Search for all the filles
	containing the component name
	"""
	filenames = glob.glob("postProcessing/actuatorLineElements/0/" + regionName + "." + componentName + "*")
	# Important to return filenames in sorted order:
	# file0, file1, ... file9, file10, ...
	return natsorted(filenames);


def positionList(regionName, componentName, ref_time):
	""" 
	Loop over all files containing
	component name and get position values
	at t = ref_time
	"""
	x, y, z = [], [], []
	filenames = getFilenames(regionName, componentName)
	index, time = getClosestTime(filenames[1], ref_time)
	for filename in filenames:
		xnew, ynew, znew = getPosition(filename, index)
		x.append(xnew)
		y.append(ynew)
		z.append(znew)

	return time, x, y, z;


def createDirectory(path):
	""" 
	If the directory does not exist,
	create it.
	"""
	if not os.path.exists(path):
		os.mkdir(path)


def writeVTK(dirName, regionName, time, x, y, z, indexes):
	""" 
	Writes VTK file with x y z coordinates. 
	Saves it to path 'dirName' 
	"""
	f = open(dirName + '/' + regionName + "_t=" + str(time) + ".vtk", "w")
	headers = ['# vtk DataFile Version 3.0', 'vtk output', 'ASCII', 'DATASET POLYDATA']
	for header in headers:
		f.write(header)
		f.write('\n')

	# Write points
	length = len(x)
	points = 'POINTS ' + str(length) + ' float \n'
	f.write(points)

	for i in range(length):
		px = str(x[i])
		py = str(y[i])
		pz = str(z[i])
		f.write(px + ' ' + py + ' ' + pz + '\n')

	# Write lines
	lines = 'LINES ' + str(length-len(indexes)-1) + ' ' + str(3*(length-len(indexes)-1)) + ' \n'
	f.write(lines)

	# Don't draw lines between different components!
	# For that, use the 'indexes' list which contains
	# the index where different components start and end
	# At that position in the list, don't write any line.

	write = 1

	for i in range(length-1):
		for num in indexes:
			if(i==num-2):
				write = 0
		if(write):
			f.write('2 ' + str(i+1) + ' ' + str(i+2) + '\n')
		write = 1

def readTimes():
	""" Open log.foamListTimes, which includes
		the time folders from OpenFOAM. Then, 
		output these times as a list.
	"""
	f = open("log.foamListTimes", "r")
	times = f.read()
	times_str = times.split("\n")
	f.close()

	del times_str[-1] # Remove last entry (is empty)
	times_float = []
	times_float.append(0) # We also want t=0
	for time in times_str:
		times_float.append(float(time))

	return times_float;


#MAIN CODE
if __name__ == "__main__":

	n = len(sys.argv)
	if n<3:
		print('Not enough arguments passed')
		components = ''
		region = ''
	else:
		region = sys.argv[1]
		components = sys.argv[2:]

	times = readTimes()
	dir_name = region + '_Geometry';
	createDirectory(dir_name)
	for time in times:
		x, y, z, componentIndex = [], [], [], []
		for component in components:
			time, x_center, y_center, z_center = positionList(region, component, time)
			# Get positions of the component
			x_comp, y_comp, z_comp = transformPoints(x_center, y_center, z_center)
			componentIndex.append(len(x_comp))
			# Combine all components in a single list
			for i in range(len(x_comp)):
				x.append(x_comp[i])
				y.append(y_comp[i])
				z.append(z_comp[i])

		# The 'componentIndex' list is used to know, in the appended lists xyz,
		# where each component geometry starts and ends 
		for i in range(1,len(componentIndex)):
			componentIndex[i] = componentIndex[i] + componentIndex[i-1] 

		writeVTK(dir_name, region, time, x, y, z, componentIndex)