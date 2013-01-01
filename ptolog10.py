#!/usr/bin/python

#Parse arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Convert p-value image to -log10(p)')
#Positional arguments
arg_parse.add_argument('pval',help='ASL image to remove outliers from. Most likely the output from mcflirt',nargs=1)
arg_parse.add_argument('out',help='Root for output image',nargs=1)
args = arg_parse.parse_args()

#Import necessary libraries
import numpy as np, nibabel as nib, sys

#Load in p-val image
try:
	pval = nib.load(args.pval[0])
except (IOError):
	print 'Error: Cannot load p-value image. Exiting...'
	sys.exit()
pvalData = np.float64(pval.get_data())

#Convert p-value image to -log10(p)
loggedData = np.ma.log10(np.ma.masked_equal(pvalData,0)) * -1

#Write out image
logged = nib.Nifti1Image(loggedData,pval.get_affine())
logged.to_filename(args.out[0] + '.nii.gz')
