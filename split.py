#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# Tyler Blazey, Summer 2011

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

#Split image into odd and even splices
even_data = input_data[:,:,0::2]
odd_data = input_data[:,:,1::2]

#Get header and setup new qform adn sform matrices
input_header = input.get_header()
old_vox = input_header.structarr['pixdim']
qform = input_header.get_qform()
sform = input_header.get_sform()
qform[2,2] = old_vox[3]*2.0
sform[2,2] = old_vox[3]*2.0

#Create header for output
out_header = nib.Nifti1Header()
out_header.set_qform(qform,code=1)
out_header.set_sform(sform,code=1)

#Write out images
even = nib.Nifti1Image(even_data,None,header=out_header)
even.to_filename(args.outroot[0] + '_1.nii.gz')
odd = nib.Nifti1Image(odd_data,None,header=out_header)
odd.to_filename(args.outroot[0] + '_2.nii.gz')


