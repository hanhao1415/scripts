#!/usr/bin/python

"""
cbf_calc.py: Python script to calculate CBF from perfusion images.

Uses the following equation from Wang et. al 2003 JMRI:

f = [ 6000 * 1000 * DeltaM * lambda ] / [ 2 * M0 * TI1 * exp(-TI2 / T1a) * qTI ]

Where:

DeltaM    = Perfusion between tag and control -> calculated from data
lambda    = blood/tissue water partition coefficient -> set constant
M0        = equilibrium brain/tissue magnetization -> acquired voxelwise
TI1       = inversion time one, labeltime -> sequence dependent
TI2       = inversion time two -> sequence dependent, varies by slice
          = T11 + w
          = TI1 + w + ( slicetime * slice number )
          = TI2 + (slicetime * slice number )
w         = delay time between saturation and excitation pulses -> sequence dependent
          = TI2 - TI1
slicetime = time to acquire each slice -> sequence dependent
          = (minTR - TI2) / # of slices
minTR     = minimum TR -> sequence dependent
qTI       = correction for difference between blood and tissue T1 and venous outflow -> set constant

Note: qTI was not used in the Wang 2003 paper. See Warmuth et. al 2003 Radiology.

Uses the filtering method from Tan et. al 2009 JMRI

Tyler Blazey, Summer 2011.

"""

#Import system modules
import sys
import argparse

#Import external modules
import numpy as np
import nibabel as nib
np.seterr(all='ignore')

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Calculate CBF from PASL perfusion images')
#Positional arguments
arg_parse.add_argument('perf',help='Path to nifti perfusion image',nargs=1)
arg_parse.add_argument('m0',help='Path to nifti m0 image',nargs=1)
arg_parse.add_argument('mask',help='Path to binary nifti brain mask',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-minTR',help='Minimum TR (ms). Default is 2260.5',type=float,
				      default=[2260.5],nargs=1)
arg_parse.add_argument('-TI1',help='First inversion interval (ms). Default is 700',type=float,
					  default=[700.0],nargs=1)
arg_parse.add_argument('-TI2',help='Second inversion interval (ms). Default is 1800',type=float,
					  default=[1800.0],nargs=1)
arg_parse.add_argument('-alp',help='Inversion efficiency. Default is 0.95',type=float,
					  default=[0.95],nargs=1)
arg_parse.add_argument('-lmbda',help='Blood/tissue water parition coefficient (ml/g). Default is \
					  0.9',type=float,default=[0.9],nargs=1)
arg_parse.add_argument('-qTI',help='Correction factor for difference between blood/tissue T1 and\
					  venous outflow. Default is 0.85',type=float,default=[0.85],nargs=1)
arg_parse.add_argument('-T1a',help='Longitudinal relxation time of blood (ms). Default is 1664.0',
					  type=float,default=[1664.0],nargs=1)
arg_parse.add_argument('-thresh',help='Standard deviation threshold for outliers. Default is 3',
					  type=float,default=[3.0],nargs=1)
arg_parse.add_argument('-unfilt',action='store_const',const=1,help='Output unfiltered CBF data.\
					  Filtering refers to removal of 3d volumes whose CBF exceeds the standard \
					  deviation by value specifed by -thresh.')
arg_parse.add_argument('-out4d',action='store_const',const=1,
					  help='Output a 4d cbf image for each output.')
args = arg_parse.parse_args()

#Load perfusion images
try:
	perf = nib.load(args.perf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find perf image at %s'%(args.perf[0])
	sys.exit()

#Load mask
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()

#Load m0
try:
	m0 = nib.load(args.m0[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find m0 image at %s'%(args.m0[0])
	sys.exit()

#Check to see images have same dimensions
if mask.get_shape() != m0.get_shape() or mask.get_shape() != perf.get_shape()[0:3]:
	print 'Mismatch in images dimensions. Check data.'
	sys.exit()

#Get image data
perf_data = perf.get_data()
mask_data = mask.get_data()
m0_data = m0.get_data()

#Create 3D and 4D masked arrays
mask_array_3d = ( mask_data - 1 ) * -1 #have to invert, as numpy masks out 1s and includes 0s
mask_array_4d = np.empty_like(perf_data)
mask_array_4d[:,:,:,:] = np.expand_dims(mask_array_3d,axis=3)

#Mask the data arrays
perf_masked_data = np.ma.array(perf_data,mask=mask_array_4d)
m0_masked_data = np.ma.array(m0_data,mask=mask_array_3d)


#Create a 3D matrix with a different scale value for each slice
sliceTime = (args.minTR[0] - args.TI2[0]) / m0_masked_data.shape[2]
scale_array = np.ma.empty_like(m0_masked_data)

#Here we only have ten slice times due to interleaved runs
for slice in range(m0_masked_data.shape[2]/2):
	acqTime = args.TI2[0] + (slice * sliceTime)
	scale = (6000.0 * 1000.0 * args.lmbda[0]) / ( 2.0 * np.exp( -acqTime / args.T1a[0]) \
			 * args.TI1[0] * args.alp[0] * args.qTI[0])
	scale_array[:,:,slice*2] = scale
	scale_array[:,:,slice*2+1] = scale

#If user wants, output unfiltered cbf images
if args.unfilt == 1:
	
	#Calculate unfiltered cbf
	cbf_unfilt_masked_data = np.ma.array(np.zeros_like(perf_masked_data),mask=mask_array_4d)
	for frame in range(cbf_unfilt_masked_data.shape[3]):
		cbf_unfilt_masked_data[:,:,:,frame] = cbf_unfilt_masked_data[:,:,:,frame] + \
						           	          np.nan_to_num((perf_masked_data[:,:,:,frame] * \
						           	          scale_array/m0_masked_data))
	#Write out unfiltered 3d cbf average
	cbf_unfilt_avg_data = np.ma.average(cbf_unfilt_masked_data,axis=3)
	cbf_unfilt_avg = nib.Nifti1Image(cbf_unfilt_avg_data,perf.get_affine())
	cbf_unfilt_avg.to_filename(args.outroot[0] + '_cbf_unfilt_avg.nii.gz')
	
	#Get and write out unfilted cbf standard deviation
	cbf_std_masked_data = np.std(cbf_unfilt_masked_data,axis=3,dtype=np.float64)
	cbf_unfilt_std = nib.Nifti1Image(cbf_std_masked_data,perf.get_affine())
	cbf_unfilt_std.to_filename(args.outroot[0] + '_cbf_unfilt_std.nii.gz')
	
	#Get unfiltered cbf variance and write it out
	cbf_unfilt_var_data = np.var(cbf_unfilt_masked_data,axis=3,dtype=np.float64)
	cbf_unfilt_var = nib.Nifti1Image(cbf_unfilt_var_data,perf.get_affine())
	cbf_unfilt_var.to_filename(args.outroot[0] + '_cbf_unfilt_var.nii.gz')
	
	#Write out unfiltered 4d cbf
	if args.out4d == 1:
		cbf_unfilt = nib.Nifti1Image(cbf_unfilt_masked_data,perf.get_affine())
		cbf_unfilt.to_filename(args.outroot[0] + '_cbf_unfilt.nii.gz')

#Get the std for each frame
frame_std_array = np.zeros(perf_masked_data.shape[3])
for frame in range(perf_masked_data.shape[3]):
	frame_std_array[frame] = np.ma.std(perf_masked_data[:,:,:,frame],dtype=np.float64)

#Check for stability
if (np.log(np.amax(frame_std_array)-np.amin(frame_std_array))) < 1:
	print "No filtering will be performed."
else:

	#Setup masks for frame and slice outliers
	frame_filt_mask = np.zeros_like(perf_data)
	slice_filt_mask = np.zeros_like(perf_data)
	
	#Get the avg and standard deviation across frames
	frame_avg = np.ma.average(perf_masked_data)
	frame_std = np.ma.std(perf_masked_data,dtype=np.float64)
	
	#Get the average std and std of stds across frames
	frame_std_avg = np.average(frame_std_array)
	frame_std_std = np.std(frame_std_array,dtype=np.float64)
	
	#Setup mean and standard deviation outlier thresholds for frames
	frame_out_avg = frame_avg + (2.5 * frame_std)
	frame_out_std = frame_std_avg + (1.5 * frame_std_std)
	
	slice_avg_array = np.zeros(perf_masked_data.shape[3])
	slice_std_array = np.zeros(perf_masked_data.shape[3])
	
	for slice in range(perf_masked_data.shape[2]):
		slice_avg = np.abs(np.ma.average(perf_masked_data[:,:,slice,:]))
 		slice_std = np.ma.std(perf_masked_data[:,:,slice,:],dtype=np.float64)
 		 
 		for frame in range(perf_masked_data.shape[3]):
 		 		
 		 		#Check for frames outliers
 		 		if slice == 0:
					if np.abs(np.ma.average(perf_masked_data[:,:,:,frame])) > frame_out_avg or frame_std_array[frame] > frame_out_std:
						frame_filt_mask[:,:,:,frame] = 1
				
				#Get the average and std for each slice with the frame
				slice_avg_array[frame] = np.abs(np.ma.average(perf_masked_data[:,:,slice,frame]))
				slice_std_array[frame] = np.ma.std(perf_masked_data[:,:,slice,frame],dtype=np.float64)
				
		#Get the average std and the std of std for the slice
		slice_std_avg = np.average(slice_std_array)
 		slice_std_std = np.std(slice_std_array,dtype=np.float64)

		#Determine the outliers thresholds for the slices
 		slice_out_avg = slice_avg + ( 2.5 * slice_std )
 		slice_out_std  = slice_std_avg + ( 1.5 * slice_std_std )
 		
 		for frame in range(perf_masked_data.shape[3]):
 			if slice_avg_array[frame] > slice_out_avg or slice_std_array[frame] > slice_out_std:
 				slice_filt_mask[:,:,slice,frame] = 1
 	
 	#Combine the frame mask and the slice mask
 	comb_filt_mask = np.ma.mask_or(frame_filt_mask,slice_filt_mask)
	
	#Calculate filtered cbf
	perf_masked_data = np.ma.array(perf_data,mask=comb_filt_mask)
	cbf_filt_masked_data = np.ma.array(np.zeros_like(perf_masked_data),mask=comb_filt_mask)

	for frame in range(cbf_filt_masked_data.shape[3]):
		cbf_filt_masked_data[:,:,:,frame] = cbf_filt_masked_data[:,:,:,frame] + \
											np.nan_to_num((perf_masked_data[:,:,:,frame] * \
											scale_array/m0_masked_data))

	#Write out 3d filtered cbf average
	cbf_filt_masked_avg_data = np.ma.average(cbf_filt_masked_data,axis=3)
	cbf_filt_masked_avg = nib.Nifti1Image(cbf_filt_masked_avg_data,perf.get_affine())
	cbf_filt_masked_avg.to_filename(args.outroot[0] + '_cbf_avg.nii.gz')

	#Get filtered standard deviation and write that out
	cbf_filt_masked_std_data = np.std(cbf_filt_masked_data,axis=3,dtype=np.float64)
	cbf_filt_masked_std = nib.Nifti1Image(cbf_filt_masked_std_data,perf.get_affine())
	cbf_filt_masked_std.to_filename(args.outroot[0] + '_cbf_std.nii.gz')

	#Get, and then write out, cbf variance
	cbf_filt_masked_var_data = np.var(cbf_filt_masked_data,axis=3,dtype=np.float64)
	cbf_filt_masked_var = nib.Nifti1Image(cbf_filt_masked_var_data,perf.get_affine())
	cbf_filt_masked_var.to_filename(args.outroot[0] + '_cbf_var.nii.gz')

	#If user wants, output filtered 4D cbf image
	if args.out4d == 1:
		cbf_filt_masked = nib.Nifti1Image(cbf_filt_masked_data,perf.get_affine())
		cbf_filt_masked.to_filename(args.outroot[0] + '_cbf.nii.gz')




