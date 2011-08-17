#!/bin/python

#import system modules
import argparse
import sys

#import external modules
import numpy as np
import nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Split a 3D image into two seperate images')
arg_parse.add_argument('input',help='Input image to be split',nargs=1)
arg_parse.add_argument('outroot',help='Root for splittled volumes',nargs=1)
args = arg_parse.parse_args()

#Load Images
try:
	input = nib.load(args.input[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading input image.'
	sys.exit()
input_data = input.get_data()

#If deterimant in less than zero, flip
if np.linalg.det(input.get_affine()) < 0:
	input_data = np.flipud(input_data)

#Split image into odd and even splices
even_data = input_data[:,:,0::2]
odd_data = input_data[:,:,1::2]

#Get reference voxel size
input_header = input.get_header()
vox = input_header.structarr['pixdim']

#Create header for output
out_header = nib.Nifti1Header()
out_header.set_qform(np.array([[vox[1],0,0,vox[1]],[0,vox[2],0,vox[2]],[0,0,vox[3]*2,
					 vox[3]*2],[0,0,0,1]]),code=0)
out_header.set_sform(np.array([[vox[1],0,0,vox[1]],[0,vox[2],0,vox[2]],[0,0,vox[3]*2,
					 vox[3]*2],[vox[1],vox[2],vox[3]*2,1]]),code=1)

#Write out images
even = nib.Nifti1Image(even_data,None,header=out_header)
even.to_filename(args.outroot[0] + '_1.nii.gz')
odd = nib.Nifti1Image(odd_data,None,header=out_header)
odd.to_filename(args.outroot[0] + '_2.nii.gz')


