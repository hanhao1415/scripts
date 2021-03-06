#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""

cbf_calc.py: Python script to calculate CBF from perfusion images.

Uses the following equation from Wang et. al 2003 JMRI:

f = [ 6000 * 1000 * DeltaM * lambda ] / [ 2 * alp* M0 * TI1 * exp(-TI2 / T1a) * qTI ]

f = [DeltaM/M0] * [ [6000*1000*lambda] / [2*alp*TI1*exp(-TI2/T1a)*qTI] ]

Where:

DeltaM    = Perfusion between tag and control -> calculated from data
lambda    = blood/tissue water partition coefficient -> set constant
M0        = equilibrium brain/tissue magnetization -> acquired voxelwise
TI1       = inversion time one, labeltime -> sequence dependent
TI2       = inversion time two -> sequence dependent, varies by slice
          = T11 + w
Acq_Time  = TI1 + w + ( slicetime * slice number )
          = TI2 + (slicetime * slice number )
w         = delay time between saturation and excitation pulses -> sequence dependent
          = TI2 - TI1
slicetime = time to acquire each slice -> sequence dependent
          = (minTR - TI2) / # of slices
minTR     = minimum TR -> sequence dependent
qTI       = correction for difference between blood and tissue T1 and venous outflow -> set constant
alp	      = Labeling efficiency -> Dependent on ASL type.

Note: qTI was not used in the Wang 2003 paper. See Warmuth et. al 2003 Radiology.

Uses the filtering method from Tan et. al 2009 JMRI

Tyler Blazey, Summer 2011.

"""

#Parse arguments
import argparse
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
arg_parse.add_argument('-mthresh',help='Standard deviation threshold for mean outliers. Default is\
					  2.5',type=float,default=[2.5],nargs=1)
arg_parse.add_argument('-sthresh',help='Standard deviation threshold for standard deviation \
					  outliers. Default is 1.5',type=float,default=[1.5],nargs=1)
arg_parse.add_argument('-unfilt',action='store_const',const=1,help='Output unfiltered CBF data.')
arg_parse.add_argument('-out4d',action='store_const',const=1,
					  help='Output a 4d cbf image for each output.')
arg_parse.add_argument('-col',action='store_const',const=1,help='Account for collated data with \
					   when accounting for slice timing.')
args = arg_parse.parse_args()

#Import system modules
import sys, numpy as np, nibabel as nib

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
try:
	m0 = nib.load(args.m0[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find m0 image at %s'%(args.m0[0])
	sys.exit()

#Check to see images have same dimensions
if mask.get_shape() != m0.get_shape() or mask.get_shape() != perf.get_shape()[0:3]:
	print 'Mismatch in image dimensions. Check data.'
	sys.exit()

#Get image data
perf_data = np.float64(perf.get_data())
mask_data = np.float64(mask.get_data())
m0_data = np.float64(m0.get_data())

#Create 3D and 4D masked arrays
mask_array_3d = ( mask_data - 1 ) * -1 #have to invert, as numpy masks out 1s and includes 0s
mask_array_4d = np.empty_like(perf_data)
mask_array_4d = np.repeat(np.expand_dims(mask_array_3d,axis=3),perf_data.shape[3],axis=3)

#Mask the data arrays
perf_masked_data = np.ma.array(perf_data,mask=mask_array_4d)
m0_masked_data = np.ma.array(m0_data,mask=mask_array_3d)


#Get number of unique slice times
if args.col == 1:
	uniq_st = m0_masked_data.shape[2] / 2
else:
	uniq_st = m0_masked_data.shape[2]

#Create a 3D matrix with a different scale value for each slice

sliceTime = (args.minTR[0] - args.TI2[0]) / ( uniq_st )
scale_array = np.ma.empty_like(m0_masked_data)
for slice in range(uniq_st):
	acqTime = args.TI2[0] + (slice * sliceTime)
	scale = (6000.0 * 1000.0 * args.lmbda[0]) / ( 2.0 * np.exp( (-1 * acqTime) / args.T1a[0]) \
			 * args.TI1[0] * args.alp[0] * args.qTI[0])	
	if args.col ==1:
		scale_array[:,:,slice*2] = scale
		scale_array[:,:,slice*2+1] = scale
	else:
		scale_array[:,:,slice] = scale

#If user wants, output unfiltered cbf images
if args.unfilt == 1:

	#Calculate unfiltered cbf
	cbf_unfilt_masked_data = np.ma.array(np.zeros_like(perf_masked_data),mask=mask_array_4d)
	for frame in range(cbf_unfilt_masked_data.shape[3]):
		cbf_unfilt_masked_data[:,:,:,frame] = np.ma.add(cbf_unfilt_masked_data[:,:,:,frame],
											  np.ma.multiply(perf_masked_data[:,:,:,frame],
											  np.ma.divide(scale_array,m0_masked_data)))
	
	#Write out unfiltered 3d cbf average
	cbf_unfilt_avg_data = np.ma.mean(cbf_unfilt_masked_data,axis=3,dtype=np.float64)
	cbf_unfilt_avg = nib.Nifti1Image(cbf_unfilt_avg_data,perf.get_affine())
	cbf_unfilt_avg.to_filename(args.outroot[0] + '_cbf_unfilt_avg.nii.gz')
	
	#Get and write out unfilted cbf standard deviation
	cbf_std_masked_data = np.ma.std(cbf_unfilt_masked_data,axis=3,dtype=np.float64)
	cbf_unfilt_std = nib.Nifti1Image(cbf_std_masked_data,perf.get_affine())
	cbf_unfilt_std.to_filename(args.outroot[0] + '_cbf_unfilt_std.nii.gz')
	
	#Get unfiltered cbf variance and write it out
	cbf_unfilt_var_data = np.ma.var(cbf_unfilt_masked_data,axis=3,dtype=np.float64)
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
	cbf_filt_masked_data = np.ma.array(np.zeros_like(perf_masked_data),mask=mask_array_4d)
else:
	
	#Setup empty masks for frame and slice outliers
	frame_filt_mask = np.ma.empty_like(perf_data)
	slice_filt_mask = np.ma.empty_like(perf_data)
	frame_outlier_count = 0
	slice_outlier_count = 0
	
	#Get the avg and standard deviation across frames
	frame_avg = np.ma.mean(perf_masked_data,dtype=np.float64)
	frame_std = np.ma.std(perf_masked_data,dtype=np.float64)

	#Get the average frame std and std of stds
	frame_std_avg = np.ma.mean(frame_std_array,dtype=np.float64)
	frame_std_std = np.ma.std(frame_std_array,dtype=np.float64)
	
	#Setup mean and standard deviation outlier thresholds for frames
	frame_out_avg = frame_avg + (args.mthresh[0] * frame_std)
	frame_out_std = frame_std_avg + (args.sthresh[0] * frame_std_std)

	#Setup arrays for slice means and standard deviations
	slice_avg_array = np.ma.zeros(perf_masked_data.shape[3])
	slice_std_array = np.ma.zeros(perf_masked_data.shape[3])
	
	for slice in range(perf_masked_data.shape[2]):
		#Get the mean and standard deviation for slice across time
		slice_avg = np.ma.abs(np.ma.mean(perf_masked_data[:,:,slice,:],dtype=np.float64))
 		slice_std = np.ma.std(perf_masked_data[:,:,slice,:],dtype=np.float64)
 		
 		for frame in range(perf_masked_data.shape[3]):
 		 		#Check for frames outliers
 		 		if slice == 0:
					if np.ma.abs(np.ma.mean(perf_masked_data[:,:,:,frame],dtype=np.float64)) > frame_out_avg \
							or frame_std_array[frame] > frame_out_std:
						frame_filt_mask[:,:,:,frame] = 1
						frame_outlier_count += 1
	
				#Get the average and std for slice within frame
				slice_avg_array[frame] = np.ma.abs(np.ma.mean(perf_masked_data[:,:,slice,frame],dtype=np.float64))
				slice_std_array[frame] = np.ma.std(perf_masked_data[:,:,slice,frame],dtype=np.float64)

		#Get the average std and the std of std for the slice
		slice_std_avg = np.ma.mean(slice_std_array)
 		slice_std_std = np.ma.std(slice_std_array,dtype=np.float64)

		#Determine the outliers thresholds for the slices
 		slice_out_avg = slice_avg + ( args.mthresh[0] * slice_std )
 		slice_out_std  = slice_std_avg + ( args.sthresh[0] * slice_std_std )

 		#Check for slice outliers
 		for frame in range(perf_masked_data.shape[3]):
 			if slice_avg_array[frame] > slice_out_avg or slice_std_array[frame] > slice_out_std:
 				slice_filt_mask[:,:,slice,frame] = 1
 				slice_outlier_count += 1
 		
 	#Combine all the masks
 	comb_filt_mask = np.ma.mask_or(frame_filt_mask,slice_filt_mask)
	comb_filt_mask = np.ma.mask_or(comb_filt_mask,np.ma.make_mask(mask_array_4d))
	brain_sum = np.sum(mask_array_4d==0); masked_sum = np.sum(comb_filt_mask==0)
	rejected = np.divide((brain_sum-masked_sum),brain_sum,dtype=np.float64) * 100
	cbf_filt_masked_data = np.ma.array(np.zeros_like(perf_masked_data),mask=comb_filt_mask)
	
	#Write out filter status
	filterFile = open(args.outroot[0] + '_filtered.txt','w')
	filterFile.write('%i frames, %i slices, %.3f percent.'%(frame_outlier_count,slice_outlier_count,rejected))
	filterFile.close()
	
#Calculate filtered cbf
for frame in range(cbf_filt_masked_data.shape[3]):
	cbf_filt_masked_data[:,:,:,frame] = np.ma.add(cbf_filt_masked_data[:,:,:,frame],
										np.ma.multiply(perf_masked_data[:,:,:,frame],
										np.ma.divide(scale_array,m0_masked_data)))

#Write out 3d filtered cbf average
cbf_filt_masked_avg_data = np.ma.mean(cbf_filt_masked_data,axis=3,dtype=np.float64)
cbf_filt_masked_avg = nib.Nifti1Image(cbf_filt_masked_avg_data,perf.get_affine())
cbf_filt_masked_avg.to_filename(args.outroot[0] + '_cbf_avg.nii.gz')

#Get filtered standard deviation and write that out
cbf_filt_masked_std_data = np.ma.std(cbf_filt_masked_data,axis=3,dtype=np.float64)
cbf_filt_masked_std = nib.Nifti1Image(cbf_filt_masked_std_data,perf.get_affine())
cbf_filt_masked_std.to_filename(args.outroot[0] + '_cbf_std.nii.gz')

#Get, and then write out, filtered cbf variance
cbf_filt_masked_var_data = np.ma.var(cbf_filt_masked_data,axis=3,dtype=np.float64)
cbf_filt_masked_var = nib.Nifti1Image(cbf_filt_masked_var_data,perf.get_affine())
cbf_filt_masked_var.to_filename(args.outroot[0] + '_cbf_var.nii.gz')

#If user wants, output filtered 4D cbf image
if args.out4d == 1:
	cbf_filt_masked = nib.Nifti1Image(cbf_filt_masked_data,perf.get_affine())
	cbf_filt_masked.to_filename(args.outroot[0] + '_cbf.nii.gz')




