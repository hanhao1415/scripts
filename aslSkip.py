#!/usr/bin/python

#Parse arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Use FSL .par files to determine frames to skip')
#Positional arguments
arg_parse.add_argument('asl',help='ASL image to remove outliers from. Most likely the output from mcflirt',nargs=1)
arg_parse.add_argument('par',help='Motion parameters estimated from mcflirt.',nargs=1)
arg_parse.add_argument('out',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-trans',help='Translation threshold. Default is 2.0 mm.',nargs=1,default=[2.0],type=float)
arg_parse.add_argument('-rot',help='Rotation threshold. Default is 2.0 degrees.',nargs=1,default=[2.0],type=float)
arg_parse.add_argument('-dtrans',help='Relative translation threshold. Default is 0.8 mm.',nargs=1,default=[0.8],type=float)
arg_parse.add_argument('-drot',help='Relative rotation threshold. Default is 0.8 degrees.',nargs=1,default=[0.8],type=float)
arg_parse.add_argument('-m0',help='Frame number for m0. Starts counting from 0.',nargs=1,default=[0])
arg_parse.add_argument('-plot',help='Output motion plots',action='store_const',const=1)
args = arg_parse.parse_args()

#Import libraries
import sys, numpy as np, nibabel as nib, matplotlib.pyplot as plt

#Read in motion parameters
try:
	motion_params = np.float64(np.loadtxt(args.par[0]))
	dims = motion_params.shape
except (IOError):
	print 'Error: Cannot load motion parameters at %s'%(args.par[0])
	sys.exit()

#Load in asl image	
try:
	asl = nib.load(args.asl[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find asl image at %s'%(args.asl[0])
	sys.exit()
	
#Check for frame number consistency
if asl.get_shape()[3] != motion_params.shape[0]:
	print 'Mismatch in number of frames between %s and %s. Check data.'%(args.par[0],args.asl[0])
	sys.exit()

#Get image data
asl_data = np.float64(asl.get_data())

#Setup dat: Convert radians to degrees and remove m0
dat = np.delete(np.concatenate((motion_params[:,0:3]*(45/np.arctan(1)),motion_params[:,3:6]),axis=1),args.m0,axis=0)

#Calculate total translations and rotations
rotSum = np.sqrt(np.power(dat[:,0],2)+np.power(dat[:,1],2)+np.power(dat[:,2],2))
transSum = np.sqrt(np.power(dat[:,3],2)+np.power(dat[:,4],2)+np.power(dat[:,5],2))

#Make outlier masks
rotOut = np.greater_equal(rotSum[:],args.rot[0])
transOut = transSum >= args.trans[0]
outMask = np.logical_or(rotOut,transOut)

#Calculate tag/control differences
diffDat = dat[::2] - dat[1::2]
diffRotSum = np.sqrt(np.power(diffDat[:,0],2)+np.power(diffDat[:,1],2)+np.power(diffDat[:,2],2));
diffTransSum = np.sqrt(np.power(diffDat[:,3],2)+np.power(diffDat[:,4],2)+np.power(diffDat[:,5],2))

#Make difference masks
diffRotOut = diffRotSum >= args.drot[0]
diffTransOut = diffTransSum >= args.dtrans[0]
diffMask = np.logical_or(diffRotOut,diffTransOut)
diffOutMask = np.empty_like(outMask); diffOutMask[::2] = diffMask; diffOutMask[1::2] = diffMask

#Combine masks
skipMask = np.logical_or(outMask,diffOutMask)

#Make evens and odds match
evenMask = skipMask[::2]; oddMask = skipMask[1::2]
combinedMask = np.logical_or(evenMask,oddMask)
skipMask[::2] = combinedMask; skipMask[1::2] = combinedMask

#Add back in m0
skipMask = np.insert(skipMask,1,0)

#Write out skip text file
np.savetxt(args.out[0]+'_skip.txt',skipMask,fmt='%i')

#Remove skipped frames from timeseries
frameAdjust = 0
for frame in range(skipMask.shape[0]):
	if skipMask[frame] == True:
		asl_data = np.delete(asl_data,frame-frameAdjust,axis=3)
		frameAdjust += 1
		
#Write out scrubbed data
asl_scrubbed = nib.Nifti1Image(asl_data,asl.get_affine())
asl_scrubbed.to_filename(args.out[0] + '_scrubbed.nii.gz')

#Make plots if user wanted
if args.plot == 1:
	#Translations Plot
	fig1 = plt.figure(1)
	trans = fig1.add_subplot(1,1,1)
	trans.plot(np.arange(1,transSum.shape[0]+1),transSum,lw=1.5,color="blue",marker="o",markersize=5)
	trans.set_xlim([0,transSum.shape[0]])
	trans.set_xlabel('Frame',fontsize=14)
	trans.set_ylabel('Translations (mm)',fontsize=14)
	trans.set_title('MCFLIRT Estimated Translations',fontsize=14)
	for frame in range(skipMask.shape[0]):
		if skipMask[frame] == True:
			trans.axvspan(frame,frame,color="red",alpha=0.5,lw=3.0)
	fig1.savefig(args.out[0] + '_trans.png')
	fig1.clf()

	#Rotation plot
	fig2 = plt.figure(2)
	rot = fig2.add_subplot(1,1,1)
	rot.plot(np.arange(1,rotSum.shape[0]+1),rotSum,lw=1.5,color="blue",marker="o",markersize=5)
	rot.set_xlim([0,rotSum.shape[0]])
	rot.set_xlabel('Frame',fontsize=14)
	rot.set_ylabel('Rotations (degrees)',fontsize=14)
	rot.set_title('MCFLIRT Estimated Rotations',fontsize=14)
	for frame in range(skipMask.shape[0]):
		if skipMask[frame] == True:
			rot.axvspan(frame,frame,color="red",alpha=0.5,lw=3.0)
	fig2.savefig(args.out[0] + '_rot.png')
	fig2.clf()







