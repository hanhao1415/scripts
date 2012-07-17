#!/bin/python

#Load in libraries
import nibabel as nib, numpy as np, argparse, sys

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Make signed negative log10 of the p-value')
#Positional arguments
arg_parse.add_argument('p',help='P-value image',nargs=1)
arg_parse.add_argument('beta',help='Signed beta image',nargs=1)
arg_parse.add_argument('outroot',help='Outroot for iamge',nargs=1)
args=arg_parse.parse_args()

#Load images
try:
        p = nib.load(args.p[0])
except (IOError,nib.spatialimages.ImageFileError):
        print 'Cannot find p-value image at %s'%(args.p[0])
        sys.exit()
pData = np.ma.masked_equal(np.float64(p.get_data()),0)

try:
	beta = nib.load(args.beta[0])
except (IOError,nib.spatialimages.ImageFileError):
        print 'Cannot find beta image at %s'%(args.beta[0])
        sys.exit()
betaData = np.float64(beta.get_data())

#Make sign
signedData = -1 * np.ma.log10(pData) * np.sign(beta)

#Output image
signed = nib.Nifti1Image(signedData,p.get_affine())
signed.to_filename(args.outroot[0] + '_logp.nii.gz')
