#!/usr/bin/python

"""

cbf_calc_3d.py: Python script to calculate CBF from average perfusion image.

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
perf_data = perf.get_data()
mask_data = mask.get_data()
m0_data = m0.get_data()

#Create mask array
mask_array = ( mask_data - 1 ) * -1 #have to invert, as numpy masks out 1s and includes 0s

#Mask the data arrays
perf_masked_data = np.ma.array(perf_data,mask=mask_array)
m0_masked_data = np.ma.array(m0_data,mask=mask_array)

#Create a 3D matrix with a different scale value for each slice
#Note, that we will assume that there are ten slices for calculation due to interleaved runs
sliceTime = (args.minTR[0] - args.TI2[0]) / ( m0_masked_data.shape[2] / 2 )
scale_array = np.ma.empty_like(m0_masked_data)
for slice in range(m0_masked_data.shape[2]/2):
	acqTime = args.TI2[0] + (slice * sliceTime)
	scale = (6000.0 * 1000.0 * args.lmbda[0]) / ( 2.0 * np.exp( (-1 * acqTime) / args.T1a[0]) \
			 * args.TI1[0] * args.alp[0] * args.qTI[0])
	scale_array[:,:,slice*2] = scale
	scale_array[:,:,slice*2+1] = scale


#Calculate cbf
cbf_masked_data = np.nan_to_num((perf_masked_data[:,:,:] * \
						 scale_array/m0_masked_data))

#Write out cbf image
cbf_masked = nib.Nifti1Image(cbf_masked_data,perf.get_affine())
cbf_masked.to_filename(args.outroot[0] + '_cbf.nii.gz')


