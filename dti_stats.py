#!/usr/bin/python

#Import system modules
import sys
import argparse

#Import external modules
import numpy as np
import nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Get DTI stats for each JHU ROI')
#Positional arguments
arg_parse.add_argument('dti',help='4D image with DTI measure at each timepoint',nargs=1)
arg_parse.add_argument('atlas',help='Atlas image where each ROI has a different intensity',nargs=1)
arg_parse.add_argument('stat',help='Name for outputed stats file',nargs=1)
args = arg_parse.parse_args()

#Load images
try:
	dti = nib.load(args.dti[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find dti image at %s'%(args.dti[0])
	sys.exit()
try:
	atlas = nib.load(args.atlas[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find atlas image at %s'%(args.atlas[0])
	sys.exit()

#Check to see images have same dimensions
if dti.get_shape()[0:3] != atlas.get_shape():
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
dtiData = np.ma.masked_equal(np.float64(dti.get_data()),0)
atlasData = atlas.get_data()

#Find number of ROIs
rois = np.unique(atlasData); rois = rois[rois!=0]

#Make empty matrix for storing results
statMatrix = np.empty((dti.shape[3],rois.shape[0]))

#Fill in stats for each ROI
for roi in range(rois.shape[0]):
	atlasMask = atlasData==rois[roi]
	statMatrix[:,roi] = np.mean(dtiData[atlasMask,:],axis=0)

#Write out result
np.savetxt(args.stat[0],statMatrix,fmt='% 2.10f')


