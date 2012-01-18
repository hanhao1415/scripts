#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# Tyler Blazey, Summer 2011

#import system modules
import sys
import argparse

#import external modules
import numpy as np
import nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Collate two seperate 4D images')
arg_parse.add_argument('img1',help='Input image 1',nargs=1)
arg_parse.add_argument('img2',help='Input image 2',nargs=1)
arg_parse.add_argument('out',help='Name for for collated image',nargs=1)
args = arg_parse.parse_args()

#Load Images
try:
	img1 = nib.load(args.img1[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading input image.'
	sys.exit()
try:
	img2 = nib.load(args.img2[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading input image.'
	sys.exit()

#Check to see input images have same dimensions
if img1.get_shape() != img2.get_shape():
	print 'Mismatch in images dimensions. Check data.'
	sys.exit()

#Load image data
img1_data = img1.get_data()
img2_data = img2.get_data()

#Create an empty image
collated_data = np.empty((img1_data.shape[0],img1_data.shape[1],img1_data.shape[2]*2,
						img1_data.shape[3]))

#Fill empty image with input data
collated_data[:,:,0::2,:] = img1_data
collated_data[:,:,1::2,:] = img2_data

#Call header and setup new qform and sform matrices
img1_header = img1.get_header()
old_vox = img1_header.structarr['pixdim']
qform = img1_header.get_qform()
sform = img1_header.get_sform()
qform[2,2] = old_vox[3]/2.0
sform[2,2] = old_vox[3]/2.0

#Create header for output
out_header = nib.Nifti1Header()
out_header.set_qform(qform,code=1)
out_header.set_sform(sform,code=1)

#Write out image
collated = nib.Nifti1Image(collated_data,None,header=out_header)
collated.to_filename(args.out[0] + '.nii.gz')