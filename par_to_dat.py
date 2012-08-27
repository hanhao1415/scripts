#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

#Parse arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Convert FSL .par file to 4dfp type .dat files')
#Positional arguments
arg_parse.add_argument('raw',help='Unaligned image. Most likely the input to mcflirt.',nargs=1)
arg_parse.add_argument('mcf',help='Realigned image. Most likely the output of mcflirt.',nargs=1)
arg_parse.add_argument('par',help='Motion parameters estimated from mcflirt.',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-mask',help='Masked used for calculating mean intensity.',nargs=1,default=[0])
arg_parse.add_argument('-rel',help='Output tracjectory relative to run mean. Same as mat2dat -R option',action='store_const',default=[0],const=[1])
arg_parse.add_argument('-diff',help='Save differentiated trajectory. Same as mat2dat -D option.',action='store_const',default=[0],const=[1])
arg_parse.add_argument('-skip',help='Number of frames to skip. Default is 3. Same as mat2dat -n option.',type=int,default=[3],nargs=1)
arg_parse.add_argument('-rad',help='Radius used to calculated total rms. Default is 50.',type=float,default=[50.0],nargs=1)
arg_parse.add_argument('-ref',help='Reference frame for scale calculation. Default is 4.',type=int,default=[4],nargs=1)
arg_parse.add_argument('-plot',help='Output motion plots',action='store_const',default=[0],const=[1])
args = arg_parse.parse_args()

#Import  modules
import sys, numpy as np, nibabel as nib, matplotlib.pyplot as plt 

#Function to calculate summary parameters
def param_sum(filename,params):
	
	#Get mean and stdev for each parameter
	mean = np.ma.mean(params,axis=0)
	stdev = np.ma.std(params,axis=0,ddof=1)
	#Get rotational and translational rms
	trans_sum = np.sum(np.power(stdev[1:4],2))
	trans_rms = np.sqrt(trans_sum)
	rot_rms = np.sqrt(np.sum(np.power(stdev[4:7],2)))
	#Get total rms
	rot_rad_sum = np.sum(np.power(np.ma.std(params[:,4:7]*(np.arctan(1)/45),axis=0,ddof=1),2))
	total_rms = np.sqrt(trans_sum+(np.power(args.rad[0],2)*rot_rad_sum))
	
	#Write out results 
	file.write('#Counting %i out of %i frames\n'%(params.shape[0]-args.skip[0],params.shape[0]))
	file.write('#mean  % .5f % .5f % .5f % .5f % .5f % .5f % .5f\n'%(mean[1],mean[2],mean[3],mean[4],mean[5],mean[6],mean[7]))
	file.write('#s.d.  % .5f % .5f % .5f % .5f % .5f % .5f % .5f\n'%(stdev[1],stdev[2],stdev[3],stdev[4],stdev[5],stdev[6],stdev[7]))
	file.write('#rms translation (mm) % .5f\n'%(trans_rms))
	file.write('#rms rotation (deg) % .5f\n'%(rot_rms))
	file.write('#total rms movement at radius=%.3fmm %.5f\n'%(args.rad[0],total_rms))
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
	trans = fig2.add_subplot(1,1,1)
	trans.plot(params[:,3:6])
	trans.set_xlim([0,dims[0]])
	trans.set_xlabel('Frame',fontsize=14)
	trans.set_ylabel('Rotations (degrees)',fontsize=14)
	trans.set_title(type + ': MCFLIRT Estimated Rotations',fontsize=14)
	pos = trans.get_position(); trans.set_position([pos.x0, pos.y0, pos.width * 0.95, pos.height])
	trans.legend(('X', 'Y', 'Z'),'upper right', shadow=True,title='Axis',bbox_to_anchor=(1.17,1))
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

#Setup image masking
if args.mask[0] == 0:
	image_mask_data = aligned_data <= 0
else:
	try:
		input_mask = nib.load(args.mask[0])
	except:
		print 'Error: Cannot load mask at %s'(args.mask[0])
		sys.exit()
	input_mask_data = ( np.float64(input_mask.get_data()) - 1 ) * -1 #have to invert, as numpy masks out 1s and includes 0s
	image_mask_data = np.empty_like(aligned_data)
	image_mask_data = np.repeat(np.expand_dims(input_mask_data,axis=3),aligned_data.shape[3],axis=3)

#Check inputs
if args.skip[0] >= aligned_data.shape[3]:
	print 'Error: Skip number cannot exceed or equal frame count.'
	sys.exit()
if raw.get_shape() != aligned.get_shape():
	print 'Error: Mismatch in image dimensions. Check data.'
	sys.exit()
if raw.get_shape()[3] != dims[0]:
	print 'Error: Motion parameters and image dimensions do not match'
	sys.exit()
if motion_params.shape[1] != 6:
	print 'Error: Motion parameters do not have the expected six columns'
	sys.exit()
if args.ref[0] <= args.skip[0]:
	print 'Error: Reference frame is within skip.'
	sys.exit()

#Mask images
raw_masked_data = np.ma.array(raw_data,mask=image_mask_data)
aligned_masked_data = np.ma.array(aligned_data,mask=image_mask_data)

#Calculate intensity scaling
ref_avg = np.ma.mean(raw_masked_data[:,:,:,args.ref[0]-1])
frame_avgs = np.empty((dims[0],1))
for frame in range(aligned_data.shape[3]): #Loop here until apply_over_axes properly handles masks
	frame_avgs[frame,:] = np.ma.mean(aligned_masked_data[:,:,:,frame])
scale = np.divide(ref_avg,frame_avgs)

#Get frame numbers
frames = np.arange(1,dims[0]+1).reshape((dims[0],1))

#Setup dat: add in frame count, convert radians to degrees and flip rotations and translations.
dat = np.concatenate((frames,motion_params[:,3:6],motion_params[:,0:3]*(45/np.arctan(1)),scale),axis=1)

#Generate masks, then apply to dat
mask = np.concatenate((np.greater_equal(np.arange(0,args.skip[0]),0),np.less(np.arange(args.skip[0],dat.shape[0]),0))).reshape(dat.shape[0],1)
params_mask = np.repeat(mask,dat.shape[1],axis=1)
dat_masked = np.ma.array(dat,mask=params_mask)

#Mask input images
#Write out a .dat type file
file = open(args.outroot[0] + '.dat','w')
file.write('#%s\n'%(' '.join(sys.argv)))
file.write('%-6s %-8s %7s %8s %8s %8s %8s %8s'%('#frame',' dX(mm)','dY(mm)','dZ(mm)','X(deg)','Y(deg)','Z(deg)','Scale\n'))
np.savetxt(file,dat,fmt=['%6i','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f'])
param_sum(file,dat_masked); file.close()
if args.plot[0] == 1: param_graph('dat',dat[:,1:7])

if args.diff[0] == 1:
	#Calculate framewise differences. Add back in first frame. Multiple scale by 100 and divide by mean.
	diff = np.diff(dat[:,1:8],axis=0); diff = np.insert(diff,0,0.0,axis=0); 
	diff[:,6] = diff[:,6] * 100.0 / np.ma.mean(dat_masked[:,7])
	
	#Writeout .ddat file
	ddat = np.concatenate((frames,diff),axis=1)
	ddat_masked = np.ma.array(ddat,mask=params_mask)
	file = open(args.outroot[0] + '.ddat','w'); 
	file.write('#%s\n'%(' '.join(sys.argv)))
	file.write('%6s %-9s %7s %8s %8s %8s %8s %9s'%('#Frame',' ddX(mm)','ddY(mm)','ddZ(mm)','dX(deg)','dY(deg)','dZ(deg)',' 100*Scale\n'))
	np.savetxt(file,ddat,fmt=['%6i','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f'])
	param_sum(file,ddat_masked); file.close()
	if args.plot[0] == 1: param_graph('ddat',ddat[:,1:7])

if args.rel[0] == 1:
	#Subtract mean motion paramaters
	param_mean = np.ma.mean(dat_masked[:,1:7],axis=0); param_rel = dat[:,1:7] - param_mean
	scale_mean = np.ma.mean(dat_masked[:,7])
	scale_rel = np.divide(dat[:,7],scale_mean).reshape(param_rel.shape[0],1)
	
	#Writeout .rdat file
	rdat = np.concatenate((frames,param_rel,scale_rel),axis=1)
	rdat_masked = np.ma.array(rdat,mask=params_mask)
	file = open(args.outroot[0] + '.rdat','w');
	file.write('#%s\n'%(' '.join(sys.argv)))
	file.write('%-6s %-8s %7s %8s %8s %8s %8s %8s'%('#frame',' dX(mm)','dY(mm)','dZ(mm)','X(deg)','Y(deg)','Z(deg)','Scale\n'))
	np.savetxt(file,rdat,fmt=['%6i','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f','% .5f'])
	param_sum(file,rdat_masked); file.close()
	if args.plot[0] == 1: param_graph('rdat',rdat[:,1:7])



