#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

#Import system modules
import sys, argparse

#Import external modules
import numpy as np, nibabel as nib, matplotlib.pyplot as plt 

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Convert FSL .par file to 4dfp type .dat files')
#Positional arguments
arg_parse.add_argument('raw',help='Unaligned image. Most likely the input to mcflirt.',nargs=1)
arg_parse.add_argument('mcf',help='Realigned image. Most likely the output of mcflirt.',nargs=1)
arg_parse.add_argument('par',help='Motion parameters estimated from mcflirt.',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-mask',help='Masked used for calculating mean intensity.',nargs=1,default=[0])
arg_parse.add_argument('-rel',help='Output tracjectory relative to run mean. Same as mat2dat -R option',action='store_const',default=[0],const=[1])
arg_parse.add_argument('-dis',help='Save differentiated trajectory. Same as mat2dat -D option.',action='store_const',default=[0],const=[1])
arg_parse.add_argument('-skip',help='Number of frames to skip. Default is 3. Same as mat2dat -n option.',type=int,default=[3],nargs=1)
arg_parse.add_argument('-rad',help='Radius used to calculated total rms. Default is 50.',type=float,default=[50.0],nargs=1)
arg_parse.add_argument('-ref',help='Reference frame for scale calculation. Default is 1.',type=int,default=[1],nargs=1)
arg_parse.add_argument('-plot',help='Output motion plots',action='store_const',default=[0],const=[1])
args = arg_parse.parse_args()

#Function to calculate summary parameters
def param_sum(filename,params):
	
	#Get mean and stdev for each parameter
	mean = np.ma.mean(params,axis=0)
	stdev = np.ma.std(params,axis=0,ddof=1)
	#Get rotational and translational rms
	dis_rms = np.ma.sqrt(np.ma.sum(np.ma.power(stdev[1:4],2)))
	rot_rms = np.ma.sqrt(np.ma.sum(np.ma.power(stdev[4:7],2)))
	#Get total rms
	total_rms = np.ma.sqrt(np.ma.power(rot_rms*args.rad[0],2)+np.ma.power(dis_rms,2)) #This is wrong currently.
	
	#Write out results 
	file.write('#Counting %i out of %i frames\n'%(params.shape[0]-args.skip[0],params.shape[0]))
	file.write('#Mean: % .5f % .5f % .5f % .5f % .5f % .5f % .5f\n'%(mean[1],mean[2],mean[3],mean[4],mean[5],mean[6],mean[7]))
	file.write('#Stdv: % .5f % .5f % .5f % .5f % .5f % .5f % .5f\n'%(stdev[1],stdev[2],stdev[3],stdev[4],stdev[5],stdev[6],stdev[7]))
	file.write('#Dis RMS: % .5f\n'%(dis_rms))
	file.write('#Rot RMS: % .5f\n'%(rot_rms))
	file.write('#Total RMS at radius %.3f: %.5f\n'%(args.rad[0],total_rms))
	return

#Function to output graphs
def param_graph(type,params):
	
	#Rotation Plot
	fig1 = plt.figure(1)
	rot = fig1.add_subplot(1,1,1)
	rot.plot(params[:,0:3])
	rot.set_xlim([0,dims[0]])
	rot.set_xlabel('Frame',fontsize=14)
	rot.set_ylabel('Translations (mm)',fontsize=14)
	rot.set_title(type + ': MCFLIRT Estimated Translations',fontsize=14)
	pos = rot.get_position(); rot.set_position([pos.x0, pos.y0, pos.width * 0.95, pos.height])
	rot.legend(('X', 'Y', 'Z'),'upper right', shadow=True,title='Axis',bbox_to_anchor=(1.17,1))
	fig1.savefig(args.outroot[0] + '_' + type + '_trans.png')
	fig1.clf()

	#Translation plot
	fig2 = plt.figure(2)
	dis = fig2.add_subplot(1,1,1)
	dis.plot(params[:,3:6])
	dis.set_xlim([0,dims[0]])
	dis.set_xlabel('Frame',fontsize=14)
	dis.set_ylabel('Rotations (degrees)',fontsize=14)
	dis.set_title(type + ': MCFLIRT Estimated Rotations',fontsize=14)
	pos = dis.get_position(); dis.set_position([pos.x0, pos.y0, pos.width * 0.95, pos.height])
	dis.legend(('X', 'Y', 'Z'),'upper right', shadow=True,title='Axis',bbox_to_anchor=(1.17,1))
	fig2.savefig(args.outroot[0] + '_' + type + '_rot.png')
	fig2.clf()

#Load raw image. 
try:
	raw = nib.load(args.raw[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error: Cannot load raw image at %s'%(args.raw[0])
	sys.exit()	
raw_data = np.float64(raw.get_data())
	
#Load realigned image.
try:
	aligned = nib.load(args.mcf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error: Cannot load realigned image at %s'%(args.mcf[0])
	sys.exit()
aligned_data = np.float64(aligned.get_data())

#Read in motion parameters
try:
	motion_params = np.float64(np.loadtxt(args.par[0]))
	dims = motion_params.shape
except (IOError):
	print 'Error: Cannot load motion parameters at %s'%(args.par[0])
	sys.exit()

#Setup masking
if args.mask[0] == 0:
	#Get rid of zeros
else
	#Load in mask

#Check inputs
if args.skip[0] >= aligned_data.shape[3]:
	print 'Error: Skip number cannot exceed or equal frame count.'
	sys.exit()
if raw.get_shape() != aligned.get_shape(): #add mask here
	print 'Error: Mismatch in image dimensions. Check data.'
	sys.exit()
if raw.get_shape()[3] != dims[0]:
	print 'Error: Motion parameters and image dimensions due not match'
	sys.exit()
if motion_params.shape[1] != 6:
	print 'Error: Motion parameters do not have the expected six columns'
	sys.exit()
if args.ref[0] >= args.skip[0]:
	print 'Error: Reference frame is within skip.'
	sys.exit()

#Calculate intensity scaling
ref_avg = np.mean(raw_data[:,:,:,args.ref[0]-1])
scale = np.divide(np.apply_over_axes(np.average,aligned_data,[0,1,2]),ref_avg).flatten().reshape((dims[0],1))

#Get frame numbers
frames = np.arange(1,dims[0]+1).reshape((dims[0],1))

#Setup dat: add in frame count, convert radians to degrees and flip rotations and translations.
dat = np.concatenate((frames,motion_params[:,3:6],motion_params[:,0:3]*57.2957795,scale),axis=1)

#Generate masks, then apply to dat
mask = np.concatenate((np.greater_equal(np.arange(0,args.skip[0]),0),np.less(np.arange(args.skip[0],dat.shape[0]),0))).reshape(dat.shape[0],1)
params_mask = np.repeat(mask,dat.shape[1],axis=1)
dat_masked = np.ma.array(dat,mask=params_mask)

#Mask input images
#Write out a .dat type file
file = open(args.outroot[0] + '.dat','w')
file.write('%s  %s\n'%('#',' '.join(sys.argv)))
file.write('%-6s %-8s %7s %8s %8s %8s %8s %8s'%('#frame',' dX(mm)','dY(mm)','dZ(mm)','X(deg)','Y(deg)','Z(deg)','Scale\n'))
np.savetxt(file,dat_masked,fmt=['%6i','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f'])
param_sum(file,dat_masked); file.close()
if args.plot[0] == 1: param_graph('dat',dat[:,1:7])

if args.dis[0] == 1:
	#Calculate framewise differences. Add back in first frame. Multiple scale by 100.
	diff = np.diff(dat[:,1:8],axis=0); diff = np.insert(diff,0,0.0,axis=0); diff[:,6] = diff[:,6] * 100.0
	
	#Writeout .ddat file
	ddat = np.concatenate((frames,diff),axis=1)
	ddat_masked = np.ma.array(ddat,mask=params_mask)
	file = open(args.outroot[0] + '.ddat','w'); 
	file.write('%s %s\n'%('#',' '.join(sys.argv)))
	file.write('%6s %-9s %7s %8s %9s %8s %8s %9s'%('#Frame',' ddX(mm)','ddY(mm)','ddZ(mm)','dX(deg)','dY(deg)','dZ(deg)',' 100*Scale\n'))
	np.savetxt(file,ddat,fmt=['%6i','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f'])
	param_sum(file,ddat_masked); file.close()
	if args.plot[0] == 1: param_graph('ddat',ddat[:,1:7])

if args.rel[0] == 1:
	#Subtract mean motion paramaters
	param_mean = np.mean(dat_masked[:,1:7],axis=0); param_rel = dat[:,1:7] - param_mean
	scale_mean = np.ma.mean(np.ma.array(np.apply_over_axes(np.mean,raw_data,[0,1,2]),mask=mask))
	scale_rel = np.divide(np.apply_over_axes(np.average,aligned_data,[0,1,2]),scale_mean).flatten().reshape((dims[0],1))
	
	#Writeout .rdat file
	rdat = np.concatenate((frames,param_rel,scale_rel),axis=1)
	rdat_masked = np.ma.array(rdat,mask=params_mask)
	file = open(args.outroot[0] + '.rdat','w');
	file.write('%s  %s\n'%('#',' '.join(sys.argv)))
	file.write('%-6s %-8s %7s %8s %8s %8s %8s %8s'%('#frame',' dX(mm)','dY(mm)','dZ(mm)','X(deg)','Y(deg)','Z(deg)','Scale\n'))
	np.savetxt(file,rdat,fmt=['%6i','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f'])
	param_sum(file,rdat_masked); file.close()
	if args.plot[0] == 1: param_graph('rdat',rdat[:,1:7])



