#!/usr/bin/python

#Import system modules
import sys,argparse,csv

#Import external modules
import numpy as np, nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Calculate weighted CBF mean for various partial volume fractions')
#Positional arguments
arg_parse.add_argument('cbf',help='CBF image',nargs=1)
arg_parse.add_argument('var',help='Variance image for weights',nargs=1)
arg_parse.add_argument('mask',help='CBF mask',nargs=1)
args=arg_parse.parse_args()

#Load images
try:
	cbf = nib.load(args.cbf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find cbf image at %s'%(args.cbf[0])
	sys.exit()
try:
	var = nib.load(args.var[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find var image at %s'%(args.var[0])
	sys.exit()
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()

#Check to see images have same dimensions
if cbf.get_shape() != var.get_shape() or cbf.get_shape() != mask.get_shape():
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
cbf_data = np.float64(cbf.get_data())
var_data = np.float64(var.get_data())
mask_data = (mask.get_data() - 1) * -1 

#Get stats file
cbf_masked_data = np.ma.array(cbf_data,mask=mask_data)
var_masked_data = np.ma.array(var_data,mask=mask_data)
mean = np.ma.average(cbf_masked_data,weights=np.ma.divide(1,var_masked_data))
print mean



	