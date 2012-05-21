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
arg_parse.add_argument('pvf',help='Gray matter partial volume fractions',nargs=1)
arg_parse.add_argument('mask',help='CBF mask',nargs=1)
arg_parse.add_argument('out',help='Output stats file',nargs=1)
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
	pvf = nib.load(args.pvf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find pvf image at %s'%(args.pvf[0])
	sys.exit()
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()

#Check to see images have same dimensions
if cbf.get_shape() != var.get_shape() or cbf.get_shape() != pvf.get_shape() or cbf.get_shape() != mask.get_shape():
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
cbf_data = np.float64(cbf.get_data())
var_data = np.float64(var.get_data())
pvf_data = pvf.get_data()
mask_data = (mask.get_data() - 1) * -1 

#Mask pvfmentation
pvf_masked_data = np.ma.array(pvf_data,mask=mask_data)

#Get stats file
stats = open(args.out[0],'w+')
low = [0.0,0.2,0.4,0.6,0.8,0.95]; high = [0.2,0.4,0.6,0.8,1.0,1.0]
for l,h in zip(low,high):
	pvf_mask = np.ma.logical_or(pvf_masked_data<=l,pvf_masked_data>=h)
	roi_cbf = np.ma.array(cbf_data,mask=pvf_mask)
	roi_weight = np.ma.divide(1,np.ma.array(var_data,mask=pvf_mask))
	roi_mean = np.ma.average(roi_cbf,weights=roi_weight)
	stats.write('%.4f\n'%(roi_mean))
stats.close()




	