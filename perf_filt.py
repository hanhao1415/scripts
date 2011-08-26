#!/usr/bin/python

"""

perf_filt.py: Python script to filter and average a ASL perfusion timeseries.

Uses the filtering method from Tan et. al 2009 JMRI

"""

#Import system modules
import sys
import argparse

#Import external modules
import numpy as np
import nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Filter and Average ASL perfusion images')
#Positional arguments
arg_parse.add_argument('perf',help='Path to nifti perfusion image',nargs=1)
arg_parse.add_argument('mask',help='Path to binary nifti brain mask',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-mthresh',help='Standard deviation threshold for mean outliers. Default is\
					  2.5',type=float,default=[2.5],nargs=1)
arg_parse.add_argument('-sthresh',help='Standard deviation threshold for standard deviation \
					  outliers. Default is 1.5',type=float,default=[1.5],nargs=1)
arg_parse.add_argument('-unfilt',action='store_const',const=1,help='Output unfiltered perfusion \
					  average.')
args = arg_parse.parse_args()
args = arg_parse.parse_args()

#Load images
try:
	perf = nib.load(args.perf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find perf image at %s'%(args.perf[0])
	sys.exit()
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()
	
#Check to see images have same dimensions
if mask.get_shape() != perf.get_shape()[0:3]:
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
perf_data = perf.get_data()
mask_data = mask.get_data()

#Create 3D and 4D masked arrays
mask_array_3d = ( mask_data - 1 ) * -1 #have to invert, as numpy masks out 1s and includes 0s
mask_array_4d = np.empty_like(perf_data)
mask_array_4d[:,:,:,:] = np.expand_dims(mask_array_3d,axis=3)

#Mask the perf image data
perf_masked_data = np.ma.array(perf_data,mask=mask_array_4d)

#If user wants, output unfiltered average perfusion image
if args.unfilt == 1:
	perf_masked_avg_data = np.ma.mean(perf_masked_data,axis=3)
	perf_masked_avg = nib.Nifti1Image(perf_masked_avg_data,perf.get_affine())
	perf_masked_avg.to_filename(args.outroot[0] + '_unfilt_avg.nii.gz')

#Get the std for each frame
frame_std_array = np.zeros(perf_masked_data.shape[3])
for frame in range(perf_masked_data.shape[3]):
	frame_std_array[frame] = np.ma.std(perf_masked_data[:,:,:,frame],dtype=np.float64)

#Check for stability
if (np.log(np.amax(frame_std_array)-np.amin(frame_std_array))) < 1:
	print "No filtering will be performed."
	perf_filt_masked_data = perf_masked_data
else:

	#Setup empty masks for frame and slice outliers
	frame_filt_mask = np.zeros_like(perf_data)
	slice_filt_mask = np.zeros_like(perf_data)
	
	#Get the avg and standard deviation across frames
	frame_avg = np.ma.mean(perf_masked_data)
	frame_std = np.ma.std(perf_masked_data,dtype=np.float64)
	
	#Get the average frame std and std of stds
	frame_std_avg = np.mean(frame_std_array)
	frame_std_std = np.std(frame_std_array,dtype=np.float64)
	
	#Setup mean and standard deviation outlier thresholds for frames
	frame_out_avg = frame_avg + (args.mthresh[0] * frame_std)
	frame_out_std = frame_std_avg + (args.sthresh[0] * frame_std_std)
	
	#Setup arrays for slice means and standard deviations
	slice_avg_array = np.zeros(perf_masked_data.shape[3])
	slice_std_array = np.zeros(perf_masked_data.shape[3])
	
	for slice in range(perf_masked_data.shape[2]):
		
		#Get the mean and standard deviation for slice across time
		slice_avg = np.abs(np.ma.mean(perf_masked_data[:,:,slice,:]))
 		slice_std = np.ma.std(perf_masked_data[:,:,slice,:],dtype=np.float64)
 		 
 		for frame in range(perf_masked_data.shape[3]):
 		 		
 		 		#Check for frames outliers
 		 		if slice == 0:
					if np.abs(np.ma.mean(perf_masked_data[:,:,:,frame])) > frame_out_avg \
							or frame_std_array[frame] > frame_out_std:
						frame_filt_mask[:,:,:,frame] = 1
				
				#Get the average and std for slice within frame
				slice_avg_array[frame] = np.abs(np.ma.mean(perf_masked_data[:,:,slice,frame]))
				slice_std_array[frame] = np.ma.std(perf_masked_data[:,:,slice,frame],
				                                  dtype=np.float64)
				
		#Get the average std and the std of std for the slice
		slice_std_avg = np.mean(slice_std_array)
 		slice_std_std = np.std(slice_std_array,dtype=np.float64)

		#Determine the outliers thresholds for the slices
 		slice_out_avg = slice_avg + ( args.mthresh[0] * slice_std )
 		slice_out_std  = slice_std_avg + ( args.sthresh[0] * slice_std_std )
 		
 		#Check for slice outliers
 		for frame in range(perf_masked_data.shape[3]):
 			if slice_avg_array[frame] > slice_out_avg or slice_std_array[frame] > slice_out_std:
 				slice_filt_mask[:,:,slice,frame] = 1
 	
 	#Combine all the masks
 	comb_filt_mask = np.ma.mask_or(frame_filt_mask,slice_filt_mask)
	comb_filt_mask = np.ma.mask_or(comb_filt_mask,np.ma.make_mask(mask_array_4d))
	perf_filt_masked_data = np.ma.array(perf_data,mask=comb_filt_mask)

#Get filtered average
perf_filt_masked_avg_data = np.ma.mean(perf_filt_masked_data,axis=3)
perf_filt_masked_avg = nib.Nifti1Image(perf_filt_masked_avg_data,perf.get_affine())
perf_filt_masked_avg.to_filename(args.outroot[0] + '_filt_avg.nii.gz')

	
	
	
	
	
	
	
	
	
	