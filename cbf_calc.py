#!/usr/bin/python

#Import system modules
import sys
import argparse

#Import external modules
import numpy as np
import nibabel as nib
np.seterr(all='ignore')

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Calculate CBF from ASL perfusion images')
#Positional arguments
arg_parse.add_argument('perf',help='Path to nifti perfusion image',nargs=1)
arg_parse.add_argument('m0',help='Path to nifti m0 image',nargs=1)
arg_parse.add_argument('mask',help='Path to nifti brain mask',nargs=1)
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
arg_parse.add_argument('-thresh',help='Standard deviation threshold for outliers. Default is 2',
					  type=float,default=[2.0],nargs=1)
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
perf_data = perf.get_data()

#Load mask
try:
	mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find mask image at %s'%(args.mask[0])
	sys.exit()
mask_data = mask.get_data()

#Load m0
try:
	m0 = nib.load(args.m0[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find m0 image at %s'%(args.m0[0])
	sys.exit()
m0_data = m0.get_data()

#Check to see if the image dimensions are the same
if m0_data.shape[0:2] != mask_data.shape[0:2] or m0_data.shape[0:2] != perf_data.shape[0:2]:
	print 'Image dimension mismatch. Check input data'
	sys.exit()

#Create 3D and 4D masked arrays
mask_array_3d = ( mask_data - 1 ) * -1 #have to invert, as numpy masks out 1s and includes 0s
mask_array_4d = np.empty_like(perf_data)
mask_array_4d[:,:,:,:] = np.expand_dims(mask_array_3d,axis=3)

#Create a 3D matrix with a different scale value for each slice
sliceTime = (args.minTR[0] - args.TI2[0]) / m0_data.shape[2]
scale_array = np.empty_like(m0_data)
for slice in range(m0_data.shape[2]):
    acqTime = args.TI2[0] + (slice * sliceTime)
    scale = (6000.0 * 1000.0 * args.lmbda[0]) / ( 2.0 * np.exp( -acqTime / args.T1a[0]) \
    		* args.TI1[0] * args.alp[0] * args.qTI[0])
    scale_array[:,:,slice] = scale
    
#Calculate cbf
cbf_masked = np.ma.array(np.zeros_like(perf_data),mask=mask_array_4d)
for frame in range(cbf_masked.shape[3]):
	cbf_masked[:,:,:,frame] = cbf_masked[:,:,:,frame] + \
						      np.nan_to_num((perf_data[:,:,:,frame]*scale_array/m0_data))

#If user wants, output unfiltered cbf images
if args.unfilt == 1:
	#Write out unfiltered 3d cbf average
	cbf_unfilt_avg_data = np.ma.average(cbf_masked,axis=3)
	cbf_unfilt_avg = nib.Nifti1Image(cbf_unfilt_avg_data,perf.get_affine())
	cbf_unfilt_avg.to_filename(args.outroot[0] + '_cbf_unfilt_avg.nii.gz')

	#Write out unfiltered 4d cbf
	if args.out4d == 1:
		cbf_unfilt = nib.Nifti1Image(cbf_masked,perf.get_affine())
		cbf_unfilt.to_filename(args.outroot[0] + '_cbf_unfilt.nii.gz')

#Create a temporal standard deviation image
cbf_stdev_masked = np.std(cbf_masked,axis=3)
cbf_stdev = nib.Nifti1Image(cbf_stdev_masked,perf.get_affine())
cbf_stdev.to_filename(args.outroot[0] + '_cbf_stdev.nii.gz')

#Find the average stdev within the brain
mean_stdev = np.ma.average(cbf_stdev_masked)

#Identify cbf outliers
up_stdev_thresh = mean_stdev + (mean_stdev * args.thresh[0])
low_stdev_thresh = mean_stdev - (mean_stdev * args.thresh[0])
outliers = []
for frame in range(cbf_masked.shape[3]):
	frame_avg = np.ma.average(cbf_masked[:,:,:,frame])
	if frame_avg >= up_stdev_thresh or frame_avg <= low_stdev_thresh:
		outliers.append(frame)
if len(outliers) == cbf_masked.shape[3]:
	print 'Number of outliers is equal to number of frames. Raise threshold and/or check data.'
	sys.exit()
cbf_masked = np.delete(cbf_masked,outliers,axis=3)

#Write out 3d filtered cbf average
cbf_filt_avg_data = np.ma.average(cbf_masked,axis=3)
cbf_filt_avg = nib.Nifti1Image(cbf_filt_avg_data,perf.get_affine())
cbf_filt_avg.to_filename(args.outroot[0] + '_cbf_avg.nii.gz')

#If user wants, output filtered, 4D cbf image
if args.out4d == 1:
	cbf_filt_avg = nib.Nifti1Image(cbf_masked,perf.get_affine())
	cbf_filt_avg.to_filename(args.outroot[0] + '_cbf.nii.gz')




