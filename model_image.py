#!/bin/python

#Load in libraries
import nibabel as nib, numpy as np, argparse, sys

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Output binary images where each model term was used')
#Positional arguments
arg_parse.add_argument('model',help='Image with base 2 numbers',nargs=1)
arg_parse.add_argument('outroot',help='Output image root',nargs=1)
args=arg_parse.parse_args()

#Load images
try:
        model = nib.load(args.model[0])
except (IOError,nib.spatialimages.ImageFileError):
        print 'Cannot find model image at %s'%(args.model[0])
        sys.exit()
modelData = np.ma.masked_equal(model.get_data(),0)

#Create empty images 
squareData = np.zeros_like(modelData); cubeData = np.zeros_like(modelData)
squareXData = np.zeros_like(modelData); cubeXData = np.zeros_like(modelData)

#Loop through each voxel. Should find a smarter way to do this.
for x in range(model.shape[0]):
        for y in range(model.shape[1]):
                for z in range(model.shape[2]):
                	if (np.ma.count(modelData[x,y,z]) == 1):
                		#Convert voxel value to binary
                		binVox = np.binary_repr(modelData[x,y,z])
                	
                		#Check to see if there are seven characters. If there isn't, quit fast.
                		if ( len(binVox) != 7 ):
                			print 'Unexpected value. Exiting...'
                			sys.exit()
                		
                		#Check for presence of each term
                		if ( binVox[2] == "1" ):
                			squareData[x,y,z] = 1
                		if ( binVox[3] == "1" ):
                			cubeData[x,y,z] = 1
                		if ( binVox[5] == "1" ):
                			squareXData[x,y,z] = 1
                		if ( binVox[6] == "1" ):
                			cubeXData[x,y,z] = 1

#Write out each image
square = nib.Nifti1Image(squareData,model.get_affine())
square.to_filename(args.outroot[0] + '_square.nii.gz')
cube = nib.Nifti1Image(cubeData,model.get_affine())
cube.to_filename(args.outroot[0] + '_cube.nii.gz')
squareX = nib.Nifti1Image(squareXData,model.get_affine())
squareX.to_filename(args.outroot[0] + '_squareX.nii.gz')
cubeX = nib.Nifti1Image(cubeXData,model.get_affine())
cubeX.to_filename(args.outroot[0] + '_cubeX.nii.gz')

                	
                	