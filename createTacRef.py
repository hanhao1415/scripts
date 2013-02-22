#!/usr/bin/python


#Create argument parser
import argparse
arg_parse = argparse.ArgumentParser(description='Create reference TAC file from two regional TAC files')
arg_parse.add_argument('tacOne',help='First input TAC file',nargs=1)
arg_parse.add_argument('tacTwo',help='Second input TAC file',nargs=1)
arg_parse.add_argument('tacOut',help='Outputted reference TAC file',nargs=1)
arg_parse.add_argument('-rsf',action='store_const',const=1,help='Use the RSF corrected column when making \
																 reference TAC. Default is uncorrected.')
args = arg_parse.parse_args()

#Load necessary libraries
import numpy as np, sys, re

#Load tacOne
try:
	tacOne = np.loadtxt(args.tacOne[0],skiprows=1)
	tacFile = open(args.tacOne[0],'r')
	oneVol = int(tacFile.readline().strip().split()[6])
	tacFile.close()
except (IOError):
	print 'Cannot read tacOne at %s. Exiting...'%(args.tacOne[0])
	sys.exit()

#Load in tacTwo
try:
	tacTwo = np.loadtxt(args.tacTwo[0],skiprows=1)
	tacFile = open(args.tacTwo[0],'r')
	twoVol = int(tacFile.readline().strip().split()[6])
	tacFile.close()
except (IOError):
	print 'Cannot read tacTwo at %s. Exiting...'%(args.tacTwo[0])
	sys.exit()

#Make sure the first three columns are the same between the files
if np.sum(tacOne[:,0:3]-tacTwo[:,0:3]) != 0:
	print 'The first three columns in %s and %s differ. Are these from the same pib_proc run? Exiting...'%(tacOne,tacTwo)
	sys.exit()

#Get a weighted average of the intensities between the two TAC files
if args.rsf == 1:
	refAvg = ((tacOne[:,4]*oneVol) + (tacTwo[:,4]*twoVol)) / (oneVol + twoVol)
	avgLabel = "Mean_(RSF)"
else:
	refAvg = ((tacOne[:,3]*oneVol) + (tacTwo[:,3]*twoVol)) / (oneVol + twoVol)
	avgLabel = "Mean"
	
#Create reference tac file
refTac = tacOne[:,0:4]; refTac[:,3] = refAvg

#Save the tac File
file = open(args.tacOut[0],'w')
file.write('%9s %15s %15s %15s %15s%10d\n'%('Frame_#','Start_Time_(s)','Duration_(s)',avgLabel,"NVoxels=",oneVol+twoVol))
for row in range(refTac.shape[0]):
	file.write('%9d %15f %15f %15f\n'%(refTac[row,0],refTac[row,1],refTac[row,2],refTac[row,3]))
file.close()