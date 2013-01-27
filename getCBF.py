#!/usr/bin/python

#Parse Arguments
import argparse 
arg_parse = argparse.ArgumentParser(description='Get group stats table from dian_asl')
#Positional arguments
arg_parse.add_argument('type',help='correction option',choices=['pvc', 'uc'],nargs=1)
arg_parse.add_argument('dir',help='directory containing dian_asl runs',nargs=1)
arg_parse.add_argument('sess',help='List of subjects in dir',nargs=1)
arg_parse.add_argument('out',help='name for outputted stats file',nargs=1)
#Optional arguments
arg_parse.add_argument('-skip',action='store_const',const=1,help='Skip subjects that don\'t data')
arg_parse.add_argument('-kernel',help='size of pvc kernel',nargs=1,default=[5],type=int)
args = arg_parse.parse_args()

#Import necessary modules
import os, sys, csv, numpy as np

#Load subject file
try:
	s = open(args.sess[0],'r')
	sessList = s.read().split()
	s.close()
except (IOError):
	print 'Cannot load %s. Exiting...'%(args.sess[0])
	sys.exit()
	
#Check to see if overall directory exists
if os.path.exists(args.dir[0]) is False:
	print 'Directory %s does not exist or you don\'t have access to it. Exiting...'%(args.dir[0])
	sys.exit()
	
#Make common name for stats files
if ( args.type[0] == "pvc" ):
	statSuffix = "_pvc" + str(args.kernel[0]) + "_cbf_avg_stats.txt"
else:
	statSuffix = "_cbf_avg_stats.txt"

#Make empty data array
roiData = np.empty([87,1])

#Loop through each session
for sess in sessList:
	
	#Load in file
	sessFile = os.path.join(args.dir[0],sess,'cbf',sess + statSuffix)
	try:
		sessData = np.loadtxt(sessFile,delimiter=',',dtype='string')
		roiData = np.hstack((roiData,np.append(sess,sessData[:,1]).reshape(87,1)))
	except (IOError):
		if (args.skip == 1):
			print 'Skipping %s'%(sess)
		else:
			print 'Cannot find %s in %s. Exiting...'%(sess,args.dir[0])
			sys.exit()
			
#Replace first column with the ROI names
roiData[:,0] = np.append('Session',sessData[:,0])

#Transpose so that each session is a row
roiData = roiData.T

#Write out stat file
try:
	np.savetxt(args.out[0],roiData,fmt="%s",delimiter=",")
except (IOError):
	print 'Cannot save ROI data at %s'%(args.out[0])
	sys.exit()


