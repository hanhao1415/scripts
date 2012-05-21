#!/usr/bin/python

#Import system and external modules
import sys,argpase,os,subprocess as sp,shutil,matplotlib.pyplot as plt

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Field Map Basis Function Optimization')
#Positional Arguments
arg_parse.add_argument('epi',help='EPI image to be unwarped',nargs=1)
arg_parse.add_argument('highres',help='High resolution structural image for registration.',nargs=1)
arg_parse.add_argument('lib',help='4D image of field map basis functions',nargs=1)
arg_parse.add_argument('mag',help='Structural image in space field maps',nargs=1)
arg_parse.add_argument('outdir',help='Output directory for results',nargs=1)
#Optional arguments
arg_parse.add_argument('-crange',help='Coefficient search ranges. Defaults are 100,000, 50,000 and 10,000',default=[100000,50000,10000],type=int,nargs='*')
arg_parse.add_argument('-drange',help='Dwell time search ranges. Defaults are 0.20, 0.10, and 0.05',default=[0.20,0.10,0.05],type=float,nargs='*')
arg_parse.add_argument('-dint',help='Initial dwell search value. Default is 0.98.',default=[0.98],type=float,nargs=1)
arg_parse.add_argument('-cint',help='Inital coefficient search value. Default is 0.',default=[0.0],type=float,nargs=1)
arg_parse.add_argument('-min',help='Minimum change in eta required for further search. Default is 0.05',default=[0.05],type=float,nargs=1)
arg_parse.add_argument('-iter',help='Number of iterations through each paramater. Default is 3',default=[3],type=int,nargs=1)

###################
###Preprocessing###
###################

#Make output directory
try:
	os.mkdir(args.outdir[0]); os.mkdir(args.outdir[0]os.chdir(args.outdir[0])
except OSError:
	print 'Directory %s already exists. Please change directory name.'%(args.outdir[0])
	sys.exit()

#Create log file
log = open('basis_opt.log','w')

#Copy over input images
src = [args.epi[0],args.highres[0],args.lib[0],args.mag[0]]; dst = ['epi','highres','lib','mag']
for s,d in zip(src,dst):
	try:
		shutil.copy(s,d)
	except IOError:
		print 'Cannot find %s image at %s. Exiting...'%(d,s)
		sys.exit()

