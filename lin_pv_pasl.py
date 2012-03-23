#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""

lin_pv_pasl.py: Python script to partial volume correct PASL images using a regression algorithm

Method taken from Asllani et. al Magnetic Resoance in Medicine 2008

Tyler Blazey, Summer 2011

"""

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
arg_parse.add_argument('mask',help='Path to binary nifti brain mask',nargs=1)
arg_parse.add_argument('csf_pbmap',help='Path to csf matter probability map',nargs=1)
arg_parse.add_argument('gm_pbmap',help='Path to gray matter probability map',nargs=1)
arg_parse.add_argument('wm_pbmap',help='Path to white matter probability map',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-kernel',help='Size of regression kernel. Default is 5',type=float,
				      default=[5.0],nargs=1)
arg_parse.add_argument('-nocsf',help='Assume that there is no perfusion in CSF.',
					  action='store_const',const=[1],default=[0])
args = arg_parse.parse_args()

#Load Images
try:
	perf = nib.load(args.perf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading perfusion image'
	sys.exit()
try:
	m0 = nib.load(args.m0[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading m0 image'
	sys.exit()
try: 
	brain_mask = nib.load(args.mask[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading brainmask image'
	sys.exit()
try:
	csf_pbmap = nib.load(args.csf_pbmap[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading csf probability map'
	sys.exit()
try:
	gm_pbmap = nib.load(args.gm_pbmap[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading gray probability map'
	sys.exit()
try:
	wm_pbmap = nib.load(args.wm_pbmap[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading white matter probability map'
	sys.exit()

#Check to see images have same dimensions
if brain_mask.get_shape() != m0.get_shape() or brain_mask.get_shape() != perf.get_shape()[0:3] or \
		brain_mask.get_shape() != csf_pbmap.get_shape() or \
		brain_mask.get_shape() != gm_pbmap.get_shape() or \
		brain_mask.get_shape() != wm_pbmap.get_shape():
	print 'Mismatch in images dimensions. Check data.'
	sys.exit()

#Load the image data
perf_data = perf.get_data()
m0_data = m0.get_data()
brain_mask_data = brain_mask.get_data()
csf_pbmap_data = csf_pbmap.get_data()
gm_pbmap_data = gm_pbmap.get_data()
wm_pbmap_data = wm_pbmap.get_data()

#Create masking arrays
mask_array_3d = (brain_mask_data - 1) * - 1 #have to invert, as numpy masks out 1s and includes 0s
mask_array_4d = np.empty_like(perf_data)
mask_array_4d = np.repeat(np.expand_dims(mask_array_3d,axis=3),perf_data.shape[3],axis=3)

#Setup masked arrays
m0_masked = np.ma.array(m0_data,mask=mask_array_3d)
perf_masked = np.ma.array(perf_data,mask=mask_array_4d)
csf_masked = np.ma.array(csf_pbmap_data,mask=mask_array_3d)
gm_masked = np.ma.array(gm_pbmap_data,mask=mask_array_3d)
wm_masked = np.ma.array(wm_pbmap_data,mask=mask_array_3d)
	
#Setup empy arrays
mcsf_pvc_data = np.ones_like(m0_data)
mgm_pvc_data = np.ones_like(m0_data)
mwm_pvc_data = np.ones_like(m0_data)
if args.nocsf[0] == 0:
	dcsf_pvc_data = np.ones_like(perf_data)
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
			if mask_array_3d[dim_1,dim_2,dim_3] == 1:
				if args.nocsf[0] == 0:
					[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
					mwm_pvc_data[dim_1,dim_2,dim_3],dcsf_pvc_data[dim_1,dim_2,dim_3,:],
					dgm_pvc_data[dim_1,dim_2,dim_3,:],dwm_pvc_data[dim_1,dim_2,dim_3,:]] = \
					[0,0,0,0,0,0]
				else:
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
				if args.nocsf[0] == 0:
					d_reg = np.column_stack((mcsf_kernel,mgm_kernel,mwm_kernel))
				else:
					d_reg = np.column_stack((mgm_kernel,mwm_kernel))
				
				#If there isn't at least three m_reg values greater than zero
				#set both dm and m voxels to zero
				if np.sum((m_reg>0)) < 3:
					if args.nocsf[0] == 0:
						[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
						mwm_pvc_data[dim_1,dim_2,dim_3],dcsf_pvc_data[dim_1,dim_2,dim_3,:],
						dgm_pvc_data[dim_1,dim_2,dim_3,:],dwm_pvc_data[dim_1,dim_2,dim_3,:]] = \
						[0,0,0,0,0,0]
					else:
						[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
						mwm_pvc_data[dim_1,dim_2,dim_3],dgm_pvc_data[dim_1,dim_2,dim_3,:],
						dwm_pvc_data[dim_1,dim_2,dim_3,:]] = [0,0,0,0,0]
				#Otherwise run a regression
				else:
					#Get the pseudoinverse for m0
					m_reg_inv = np.linalg.pinv(m_reg)
					
					#Get the m0 values and run a regression with them
					m0_kernel = np.ma.compressed(m0_masked[x_north:x_south,y_west:y_east,dim_3])
					[mcsf_pvc_data[dim_1,dim_2,dim_3],mgm_pvc_data[dim_1,dim_2,dim_3],
					mwm_pvc_data[dim_1,dim_2,dim_3]] = np.dot(m_reg_inv,m0_kernel)
					
					if args.nocsf[0] == 1 and np.sum((d_reg>0)) < 2:
						[dgm_pvc_data[dim_1,dim_2,dim_3,:],dwm_pvc_data[dim_1,dim_2,dim_3,:]] = \
						[0,0]
					
					else:
						#Get the pseudoinverse for perfusion
						d_reg_inv = np.linalg.pinv(d_reg)
					
						#Loop through every perf dimension
						for dim_4 in range(perf.shape[3]):
						
							#Get the perf values for each dim4 and run a regression with them
							perf_kernel = np.ma.compressed(perf_masked[x_north:x_south,
													   	   y_west:y_east,dim_3,dim_4])
						
							if args.nocsf[0] == 0:
								[dcsf_pvc_data[dim_1,dim_2,dim_3,:],
								dgm_pvc_data[dim_1,dim_2,dim_3,dim_4],
								dwm_pvc_data[dim_1,dim_2,dim_3,dim_4]] = \
								np.dot(d_reg_inv,perf_kernel)
							else:
								[dgm_pvc_data[dim_1,dim_2,dim_3,dim_4],
								dwm_pvc_data[dim_1,dim_2,dim_3,dim_4]] = \
								np.dot(d_reg_inv,perf_kernel)
	print 'Finished processing slice %s'%(dim_3+1)


#Write out some results
mcsf_pv = nib.Nifti1Image(mcsf_pvc_data,m0.get_affine())
mcsf_pv.to_filename(args.outroot[0] + '_mcsf.nii.gz')
mgm_pvc = nib.Nifti1Image(mgm_pvc_data,m0.get_affine())
mgm_pvc.to_filename(args.outroot[0] + '_mgm.nii.gz')
mwm_pvc = nib.Nifti1Image(mwm_pvc_data,m0.get_affine())
mwm_pvc.to_filename(args.outroot[0] + '_mwm.nii.gz')
if args.nocsf[0] == 0:
	dcsf_pvc = nib.Nifti1Image(dcsf_pvc_data,perf.get_affine())
	dcsf_pvc.to_filename(args.outroot[0] + '_dcsf.nii.gz')
dgm_pvc = nib.Nifti1Image(dgm_pvc_data,perf.get_affine())
dgm_pvc.to_filename(args.outroot[0] + '_dgm.nii.gz')
dwm_pvc = nib.Nifti1Image(dwm_pvc_data,perf.get_affine())
dwm_pvc.to_filename(args.outroot[0] + '_dwm.nii.gz')










