#!/data/nil-bluearc/benzinger2/Tyler/epd-7.2-2-rh5-x86_64/bin/python
#Convert between FSL style bvals and bvecs and 4dfp param files.

#Parse arguments
import argparse 
arg_parse = argparse.ArgumentParser(description='Convert FSL style bvals/bvecs to 4dfp param format')
#Positional arguments
arg_parse.add_argument('bvals',help='FSL bvals file',nargs=1)
arg_parse.add_argument('bvecs',help='FSL bvecs file',nargs=1)
arg_parse.add_argument('param',help='4dfp format output',nargs=1)
arg_parse.add_argument('-flip',help='Flip sign of y and z axes',action='store_const',const=1)
args = arg_parse.parse_args()

#Import necessary libraries
import numpy as np, sys

#Load in bvals
try:
	bvals = np.loadtxt(args.bvals[0])
except (IOError):
	print 'Cannot find bvals file at %s'%(args.bvals[0])
	sys.exit()

#Load in bvecs
try:
	bvecs = np.loadtxt(args.bvecs[0])
except (IOError):
	print 'Cannot find bvecs file at %s'%(args.bvecs[0])
	sys.exit()

#Make sure dimensions match
if bvecs.shape[1] != bvals.shape[0] or bvecs.shape[0] != 3:
	print 'Unexpected dimensions for bvals/bvecs. Check input.'
	sys.exit()

#Make param matrix
param = np.hstack((bvals.reshape(bvals.shape[0],1),bvecs.T))

#Add negative signs to y and z if asked to
if ( args.flip == 1 ):
	param[:,2] = param[:,2] * -1
	param[:,3] = param[:,3] * -1

#Write out new file
file = open(args.param[0],'w')
file.write('%i\n'%(bvals.shape[0]))
np.savetxt(args.param[0],param,fmt='%4i % 1.4f % 1.4f % 1.4f')
file.close()
