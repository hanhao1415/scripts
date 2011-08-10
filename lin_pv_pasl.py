#!/usr/bin/python

#import system modules
import sys
import argparse

#import external modules
import numpy as np
import nibabel as nib


#Parse arguments
arg_parse = argparse.ArgumentParser(description='Linear regression technique for ASL partial \
									volume correction')
#Positional arguments
arg_parse.add_argument('perf',help='Path to nifti perfusion image',nargs=1)
arg_parse.add_argument('m0',help='Path to nifti m0 image',nargs=1)
arg_parse.add_argument('mask',help='Path to nifti brain mask',nargs=1)
arg_parse.add_argument('csf_pvm',help='Path to csf matter partial volume map',nargs=1)
arg_parse.add_argument('gm_pvm',help='Path to gray matter partial volume map',nargs=1)
arg_parse.add_argument('wm_pvm',help='Path to white matter partial volume map',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-kernel',help='Size of regression kernel',type=float,
				      default=[5.0],nargs=1)
args = arg_parse.parse_args()

#Load Images
try:
	perf = nib.load(args.perf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error perf'
	sys.exit()
perf_data = perf.get_data()
try:
	m0 = nib.load(args.m0[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error m0'
	sys.exit()
m0_data = m0.get_data()
try: 
	brain_mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error brainmask'
	sys.exit()
brain_mask_data = brain_mask.get_data()
try:
	csf_pvmap = nib.load(args.csf_pvm[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error csf'
	sys.exit()
csf_pvmap_data = csf_pvmap.get_data()
try:
	gm_pvmap = nib.load(args.gm_pvm[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error gray'
	sys.exit()
gm_pvmap_data = gm_pvmap.get_data()
try:
	wm_pvmap = nib.load(args.wm_pvm[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error for wm'
	sys.exit()
wm_pvmap_data = wm_pvmap.get_data()

#Create masking arrays
mask_array = (brain_mask_data - 1) * - 1
mask_array_4d = np.empty_like(perf_data)
mask_array_4d[:,:,:,:] = np.expand_dims(mask_array,axis=3)

#Setup masked arrays
m0_masked = np.ma.array(m0_data,mask=mask_array)
perf_masked = np.ma.array(perf_data,mask=mask_array_4d)
csf_masked = np.ma.array(csf_pvmap_data,mask=mask_array)
gm_masked = np.ma.array(gm_pvmap_data,mask=mask_array)
wm_masked = np.ma.array(wm_pvmap_data,mask=mask_array)
	
#Setup empy arrays
mcsf_pvc_data = np.ones_like(m0_data)
mgm_pvc_data = np.ones_like(m0_data)
mwm_pvc_data = np.ones_like(m0_data)
dgm_pvc_data = np.ones_like(perf_data)
dwm_pvc_data = np.ones_like(perf_data)

#Setup step size
x_dim,y_dim,z_dim = m0_data.shape
step_size = (0.5 * args.kernel[0]) - 0.5

#Loop through each voxel
for dim_3 in range(perf.shape[2]):
	for dim_2 in range(perf.shape[1]):
		for dim_1 in range(perf.shape[0]):	
			
			#If voxel is outside of brain mask, set that voxel to zero for all maps
			if mask_array[dim_1,dim_2,dim_3] == 1:
				[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
				mwm_pvc_data[dim_1,dim_2,dim_3],dgm_pvc_data[dim_1,dim_2,dim_3,:],
				dwm_pvc_data[dim_1,dim_2,dim_3,:]] = [0,0,0,0,0]
			
			#Check to see if there is enough values for a regression
			else:
				
				#Setup general boundries
				x_north = dim_1 - step_size;
				x_south = dim_1 + step_size + 1;
				y_west = dim_2 - step_size;
				y_east = dim_2 + step_size + 1;
	
				#Fix step sizes if they go over matrix boundries
				if x_north < 0: x_north = 0
				if x_south > x_dim: x_south = x_dim 
				if y_west < 0: y_west = 0
				if y_east > y_dim: y_east = y_dim
		
				#Setup regression vectors and matricies with non-masked data
				mcsf_kernel = np.ma.compressed(csf_masked[x_north:x_south,y_west:y_east,dim_3])
				mgm_kernel = np.ma.compressed(gm_masked[x_north:x_south,y_west:y_east,dim_3])
				mwm_kernel = np.ma.compressed(wm_masked[x_north:x_south,y_west:y_east,dim_3])
				m_reg = np.column_stack((mcsf_kernel,mgm_kernel,mwm_kernel))
				dm_reg = np.column_stack((mgm_kernel,mwm_kernel))
				
				#If there isn't at least three m_reg values greater than zero
				#set both dm and m voxels to zero
				if np.sum((m_reg>0)) < 3:
					[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
					mwm_pvc_data[dim_1,dim_2,dim_3],dgm_pvc_data[dim_1,dim_2,dim_3,:],
					dwm_pvc_data[dim_1,dim_2,dim_3,:]] = [0,0,0,0,0]
				#Otherwise run a regression
				else:
					#Get the pseudoinverse
					m_reg_inv = np.linalg.pinv(m_reg)
					
					#Get the m0 values and run a regression with them
					m0_kernel = np.ma.compressed(m0_masked[x_north:x_south,y_west:y_east,dim_3])
					[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
					mwm_pvc_data[dim_1,dim_2,dim_3]] = np.dot(m_reg_inv,m0_kernel)
				
					#If there isn't a least three dm_reg values greater than zero, 
					#set dm voxels to zero
					if np.sum((dm_reg>0)) < 3:
						[dgm_pvc_data[dim_1,dim_2,dim_3,:],
						 dwm_pvc_data[dim_1,dim_2,dim_3,:]] = [0,0]
					#Otherwise run a regression
					else:
						#Get the pseudoinverse
						dm_reg_inv = np.linalg.pinv(dm_reg)
					
						#Loop through every perf dimension
						for dim_4 in range(perf.shape[3]):
						
							#Get the perf values for each dim4 and run a regression with them
							perf_kernel = np.ma.compressed(perf_masked[x_north:x_south,
														   y_west:y_east,dim_3,dim_4])
							[dgm_pvc_data[dim_1,dim_2,dim_3,dim_4],
							 dwm_pvc_data[dim_1,dim_2,dim_3,dim_4]] = np.dot(dm_reg_inv,perf_kernel)
	print 'Finished processing slice %s'%(dim_3+1)


#Write out some results
mgm_pvc = nib.Nifti1Image(mgm_pvc_data,m0.get_affine())
mgm_pvc.to_filename(args.outroot[0] + '_mgm.nii.gz')
mwm_pvc = nib.Nifti1Image(mwm_pvc_data,m0.get_affine())
mwm_pvc.to_filename(args.outroot[0] + '_mwm.nii.gz')
mcsf_pv = nib.Nifti1Image(mcsf_pvc_data,m0.get_affine())
mcsf_pv.to_filename(args.outroot[0] + '_mcsf.nii.gz')
dgm_pvc = nib.Nifti1Image(dgm_pvc_data,perf.get_affine())
dgm_pvc.to_filename(args.outroot[0] + '_dgm.nii.gz')
dwm_pvc = nib.Nifti1Image(dwm_pvc_data,perf.get_affine())
dwm_pvc.to_filename(args.outroot[0] + '_dwm.nii.gz')










