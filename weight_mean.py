#!/usr/bin/python

#Import system modules
import sys,argparse,csv

#Import external modules
import numpy as np, nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Calculate weighted CBF mean for each Freesurfer ROI')
#Positional arguments
arg_parse.add_argument('cbf',help='CBF image',nargs=1)
arg_parse.add_argument('var',help='Variance image for weights',nargs=1)
arg_parse.add_argument('seg',help='ROI segmentation image',nargs=1)
arg_parse.add_argument('mask',help='CBF mask',nargs=1)
arg_parse.add_argument('rois',help='Lookup table for ROI image',nargs=1)
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
	seg = nib.load(args.seg[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find seg image at %s'%(args.seg[0])
	sys.exit()
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()

#Check to see images have same dimensions
if cbf.get_shape() != var.get_shape() or cbf.get_shape() != seg.get_shape() or cbf.get_shape() != mask.get_shape():
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
cbf_data = np.float64(cbf.get_data())
var_data = np.float64(var.get_data())
seg_data = seg.get_data()
mask_data = (mask.get_data() - 1) * -1 

#Load roi file
try:
	rois = csv.reader(open(args.rois[0],'r'))
except(IOError):
	print 'Cannot open roi file at %s'%(args.rois[0])
	sys.exit()

#Mask segmentation
seg_masked_data = np.ma.array(seg_data,mask=mask_data)

#Get stats file
stats = open(args.out[0],'w+')
for roi in rois:
	print roi[1]
	roi_mask = seg_masked_data!=int(roi[0])
	roi_cbf = np.ma.array(cbf_data,mask=roi_mask)
	roi_weight = np.ma.divide(1,np.ma.array(var_data,mask=roi_mask))
	roi_mean = np.ma.average(roi_cbf,weights=roi_weight)
	stats.write('%s,%.4f\n'%(roi[1],roi_mean))
stats.close()




	