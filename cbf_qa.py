#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

#Import system modules
import sys
import argparse

#Import external modules
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Take CBF QA Screenshots')
arg_parse.add_argument('cbf',help='Path to nifti cbf image',nargs=1)
arg_parse.add_argument('out',help='Output image. Include extension',nargs=1)
arg_parse.add_argument('-ratios',help='Set of three slice ratios',nargs=4,default=[0.35,0.45,0.55,0.65])
args = arg_parse.parse_args()

#Load images
try:
	cbf = nib.load(args.cbf[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find cbf image at %s'%(args.cbf[0])
	sys.exit()
cbf_data = np.float64(cbf.get_data())

#Set plotting defaults
plt.rc('text',**{'color':'white'});
plt.rc('ytick',**{'color':'white','labelsize':'small'})

#Setup red/blue colormap
cdict = {'red': ((0.,0.,0.),(0.075,.0,0.0),(0.1499,0.,0.),(.45,.9,.9),(1,1,1)),
		 'green': ((0.,0.,0.),(0.075,.0,0.0),(0.1499,0.,0.),(.45,0,0),(1,1,1)),
		 'blue': ((0.,1,1),(0.075,.75,0.75),(0.1499,0.25,0.),(.45,0,0),(1,0,0))}
my_cmap = plt.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

#Setup slice ranges
ratios = np.hstack([args.ratios,args.ratios,args.ratios]); dims = cbf_data.shape

#Configure figure
fig = plt.figure(figsize=(20, 2), dpi=100)
gs = grd.GridSpec(1,12,width_ratios=[.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.6])

#Plot images
for i in range(0,12):
	ax = plt.subplot(gs[i])
	ax.axes.get_xaxis().set_visible(False)
	ax.axes.get_yaxis().set_visible(False)
	if i <= 3:
		slice = np.round(ratios[i]*dims[0])
		ax.imshow(np.rot90(cbf_data[:,:,slice]),cmap=my_cmap,vmin=-20,vmax=100,interpolation="bicubic")
	elif i > 3 and i <= 7:
		slice = np.round(ratios[i]*dims[1])
		ax.imshow(np.rot90(cbf_data[:,slice,:]),cmap=my_cmap,vmin=-20,vmax=100,interpolation="bicubic")
	elif i > 7 and i <= 11:
		slice = np.round(ratios[i]*dims[2])
		img = ax.imshow(np.rot90(cbf_data[slice,:,:]),cmap=my_cmap,vmin=-20,vmax=100,interpolation="bicubic")
		if i == 11:
			cbar = fig.colorbar(img,ticks=[-20,20,60,100])
fig.savefig(args.out[0],facecolor='black',bbox_inches='tight')
