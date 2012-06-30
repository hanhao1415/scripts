#!/usr/bin/python 
#Code shamelessly copied from Tom Nichol's FDR.m

#Import system modules
import sys
import argparse

#Import external modules
import numpy as np
import nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='FDR Correction for Volume Images')
#Positional arguments
arg_parse.add_argument('raw',help='Image consisiting of signed -log(p) values from a two-tailed test',nargs=1)
arg_parse.add_argument('mask',help='Image mask',nargs=1)
arg_parse.add_argument('out',help='FDR corrected image',nargs=1)
#Optional arguments
arg_parse.add_argument('-q',help='FDR threshold. Default is 0.05',type=float,default=[0.05],nargs=1)
args = arg_parse.parse_args()

#Load images
try:
	raw = nib.load(args.raw[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find raw image at %s'%(args.raw[0])
	sys.exit()
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()

#Check to see images have same dimensions
if raw.get_shape() != mask.get_shape():
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
rawData = np.float64(raw.get_data())
maskData = np.float64(mask.get_data())

#Mask data
rawMasked = np.ma.power(10,-1*np.abs(rawData[np.ma.make_mask(maskData)]))

#Sort absolute p-values smallest to largest
p = np.ma.sort(rawMasked)

#Get array length and index array
V = p.shape[0]; I = np.arange(1,V+1,dtype=np.float64)

#Calculate FDR threshold
outData = rawData; index = p[p<=(I/V)*args.q[0]]
if np.ma.sum(index) == 0:
	outData[:,:,:] = 0
else:
	thresh = -1 * np.log10(np.ma.amax(index))
	outData[np.abs(outData)<=thresh] = 0

#Output image
out = nib.Nifti1Image(outData,raw.get_affine())
out.to_filename(args.out[0] + '.nii.gz')


