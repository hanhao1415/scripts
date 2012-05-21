#!/usr/bin/python

#Import system and external modules
import sys,argparse,os,subprocess as sp,shutil,matplotlib.pyplot as plt

###########
###Usage###
###########

arg_parse = argparse.ArgumentParser(description='Field Map Basis Function Optimization')
#Positional Arguments
arg_parse.add_argument('epi',help='EPI image to be unwarped',nargs=1)
arg_parse.add_argument('highres',help='High resolution structural image for registration.',nargs=1)
arg_parse.add_argument('lib',help='4D image of field map basis functions',nargs=1)
arg_parse.add_argument('mean',help='3D mean field map image',nargs=1)
arg_parse.add_argument('mag',help='Structural image in space field maps',nargs=1)
arg_parse.add_argument('magmask',help='Binary brain mask for mag image',nargs=1)
arg_parse.add_argument('outdir',help='Output directory for results',nargs=1)
#Optional arguments
arg_parse.add_argument('-crange',help='Coefficient search ranges. Defaults are 100,000, 50,000 and 10,000',default=[100000,50000,10000],type=int,nargs='*')
arg_parse.add_argument('-drange',help='Dwell time search ranges. Defaults are 0.20, 0.10, and 0.05',default=[0.20,0.10,0.05],type=float,nargs='*')
arg_parse.add_argument('-dint',help='Initial dwell search value. Default is 0.98.',default=[0.98],type=float,nargs=1)
arg_parse.add_argument('-cint',help='Inital coefficient search value. Default is 0.',default=[0.0],type=float,nargs=1)
arg_parse.add_argument('-min',help='Minimum change in eta required for further search. Default is 0.05',default=[0.05],type=float,nargs=1)
arg_parse.add_argument('-iters',help='Number of iterations through each paramater. Default is 3',default=[3],type=int,nargs=1)
args = arg_parse.parse_args()

###################
###Preprocessing###
###################

#Make output directories
try:
	pdir = os.path.join(args.outdir[0],'preproc')
	os.mkdir(args.outdir[0]); os.mkdir(pdir);
except OSError:
	print 'Directory %s already exists. Please change directory name.'%(args.outdir[0])
	sys.exit()

#Create log file
log = open(os.path.join(args.outdir[0],'basis_opt.log'),'w')

#Convert all images to nifti format
src = [args.epi[0],args.highres[0],args.lib[0],args.mean[0],args.mag[0],args.magmask[0]]
dst = ['epi','highres','lib','mean','mag','mag_mask']
for s,d in zip(src,dst):
	try:
		convert = sp.call('mri_convert %s %s.nii'%(s,os.path.join(pdir,d)),shell=True,stdout=log,stderr=log)
	except IOError:
		print 'Cannot find %s image at %s. Exiting...'%(d,s)
		sys.exit()
os.chdir(pdir)

#Set output type to nii to appease nifti_4dfp
os.putenv('FSLOUTPUTTYPE','NIFTI')

#Get the number of basis maps in library
val = sp.Popen('fslval lib dim4',shell=True,stdout=sp.PIPE,stderr=log)
dim = val.communicate()

#Split the basis library into seperate maps
split = sp.call('fslsplit lib basis -t',shell=True,stdout=log,stderr=log)

#Register epi and field map mag image
reg = sp.call('flirt -in epi -ref mag -dof 6 -cost mutualinfo -omat epi_to_mag.mat -out epi_to_mag',shell=True,stdout=log,stderr=log)

#Convert epi and highres to 4dfp image
for image in ['epi','highres']:
	convert = sp.call('nifti_4dfp -4 %s %s'%(image,image),shell=True,stdout=log,stderr=log)

#Get an inital registration between epi and highres
modes = [1027,1027,3083]
for mode in modes:
	imgreg = sp.call('imgreg_4dfp highres none epi none epi_to_highres_t4 %i'%(mode),shell=True,stdout=log,stderr=log)
	
#for iter in range(args.iters[0]):
#	idir = os.path.join(args.outdir[0],'iter%i'%(iter)); os.mkdir(idir); os.chdir(idir)
#	params = ['dwell',
#	for 
	
