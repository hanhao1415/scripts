#!/usr/bin/python
"""

Save medial and lateral views of FreeSurfer surface.
Image pasting code modified from http://29a.ch/2009/5/14/concatenating-images-using-python

"""

#Parse some arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Combine Pysurfer Views')
#Positional arguments
arg_parse.add_argument('data',help='Statistical map in FreeSurfer surface space.',nargs=1)
arg_parse.add_argument('hemi',help='Freesurfer hemisphere',nargs=1)
arg_parse.add_argument('out',help='Name of output image.',nargs=1)
#Optional arguments
arg_parse.add_argument('-surf',help='Freesurfer surface to display map on. Default is pial.',nargs=1,default=['pial'])
arg_parse.add_argument('-minmax',help='Minumum and maximum for map. Default is 2 and 5.',nargs=2,default=[2,5],type=float)
arg_parse.add_argument('-sign',help='Sign for overlay. By default it is determined from data',nargs=1,default=[0])
arg_parse.add_argument('-negsign',help='Override sign specifed by -minmax for negative bar.',nargs=2,default=[0,0],type=float)
args = arg_parse.parse_args()

#Necessary libraries
import surfer as surf, Image, os, numpy as np, matplotlib.pyplot as plt, matplotlib.colors as col
from mayavi import mlab

#Move both colorbars
def setOverlay(overlay):
	
	#Get the data for the autumn colormap
	autumnData = plt.cm.autumn(np.arange(255))
	
	#Make a colormap with inverted data
	invertData = np.ones((255,4))
	invertData[:,0:3] = autumnData[:,0:3] * -1 + 1
	
	#Change positive colorbar
	if ( sign == 'abs' or sign == 'pos' ):
		overlay.pos_bar.scalar_bar_representation.position2 = (0.15,0.45)
		overlay.pos_bar.scalar_bar_representation.position = (0.82,0.55)
		overlay.pos_bar.scalar_bar_representation.orientation = 1
		overlay.pos_bar.lut.table = autumnData * 255
	#Change negative colorbar
	if ( sign == 'abs' or sign == 'neg' ):
		overlay.neg_bar.scalar_bar_representation.position2 = (0.15,0.45)
		overlay.neg_bar.scalar_bar_representation.position = (0.82,0.075)
		overlay.neg_bar.scalar_bar_representation.orientation = 1
		overlay.neg_bar.lut.table = invertData * 255
		if ( args.negsign[0] != 0 or args.negsign[1] != 0 ):
			overlay.neg_bar.data_range = np.array([args.negsign[0],args.negsign[1]])

#Convert float to rounded integer
def intRound(float):
	return(int(round(float)))

pics = []; outDir = os.path.dirname(args.out[0])

#Get brain image up
brain = surf.Brain('fsaverage',args.hemi[0],args.surf[0],config_opts={'size':'800'})

#Read in data
data = surf.io.read_scalar_data(args.data[0])

#Theshold the data by the minimum
data[np.abs(data)<args.minmax[0]] = 0
	
#Get minimum and maximum of data
dataMax = np.amax(data); dataMin = np.amin(data)
	
#Only show overlay if there is a need to
if ( dataMax != 0 or dataMin != 0 ):

	#Determine sign of overlay display
	if args.sign[0] != 0:
		sign = args.sign[0]
	elif dataMin < 0 and dataMax > 0:
		sign = "abs"
	elif dataMin < 0 and dataMax == 0:
		sign = "neg"
	else:
		sign = "pos"
	
	#Add overlay
	brain.add_overlay(data,min=args.minmax[0],max=args.minmax[1],sign=sign,name='multi')
	setOverlay(brain.overlays['multi'])

#Take a screenshots for each view
for view in ['m','l']:
	#Show appropriate view and adjust placement
	brain.show_view(view)
	brain.show_view({'distance': 475})
	mlab.move(0,25,0)
	
	#Take screeshot
	pic = os.path.join(outDir,'temp_%s.jpg'%(view))
	brain.save_image(pic)
	pics.append(pic)

#Read in each image
images = map(Image.open, pics)

#Determine image dimensions. Assume that all are the same
width = min(img.size[0] for img in images)
height = min(img.size[1] for img in images)

#Crop each image
cropped = []
for img in range(len(images)):
	#Crop out colorbar
	if img != 1:
		coords = (0,0,intRound(width*0.75),height)
	else:
		coords = (intRound(width*.075),0,width,height)
	crop = images[img].crop(coords)
	cropped.append(crop)

#Make new image with a width of the cropped images
cropWidth = sum(img.size[0] for img in cropped)
output = Image.new("RGB",(cropWidth,height))

#Join images along horizontal axis
startWidth = 0
for img in cropped:
	output.paste(img,(startWidth,0))
	startWidth += img.size[0]
output.save(args.out[0])

#Remove the temp files
for pic in pics:
	os.remove(pic)