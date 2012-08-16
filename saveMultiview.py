#!/usr/bin/python
"""

Save multiple for image views of FreeSurfer surface.
Image pasting code modified from http://29a.ch/2009/5/14/concatenating-images-using-python

"""

#Parse some arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Combine Pysurfer Views')
#Positional arguments
arg_parse.add_argument('map',help='Statistical map in FreeSurfer surface space.',nargs=1)
arg_parse.add_argument('hemi',help='FreeSurfer hemisphere.',nargs=1)
arg_parse.add_argument('out',help='Name of output image.',nargs=1)
#Optional arguments
arg_parse.add_argument('-surf',help='Freesurfer surface to display map on. Default is pial.',nargs=1,default=['pial'])
arg_parse.add_argument('-minmax',help='Minumum and maximum for map. Default is 2 and 5.',nargs=2,default=[2,5],type=float)
arg_parse.add_argument('-sign',help='Sign for overlay. By default it is determined form data',nargs=1,default=[0])
args = arg_parse.parse_args()

#Necessary libraries
import surfer as surf, Image, os
from mayavi import mlab
import numpy as np

#Move both colorbars
def setOverlay(overlay,sign):
	if ( sign == 'abs' or sign == 'pos' ):
		overlay.pos_bar.scalar_bar_representation.position2 = (0.15,0.45)
		overlay.pos_bar.scalar_bar_representation.position = (0.82,0.55)
		overlay.pos_bar.scalar_bar_representation.orientation = 1
		overlay.pos_bar = True
	if ( sign == 'abs' or sign == 'neg' ):
		overlay.neg_bar.scalar_bar_representation.position2 = (0.15,0.45)
		overlay.neg_bar.scalar_bar_representation.position = (0.82,0.075)
		overlay.neg_bar.scalar_bar_representation.orientation = 1

#Convert float to rounded integer
def intRound(float):
	return(int(round(float)))

#Get brain image up
brain = surf.Brain('fsaverage',args.hemi[0],args.surf[0],config_opts={'size':'800'})

#Read in data
data = surf.io.read_scalar_data(args.map[0])
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

	#Show overlay
	brain.add_overlay(data,min=args.minmax[0],max=args.minmax[1],sign=sign,name='multi')
	setOverlay(brain.overlays['multi'],sign)

#Take a screenshots for each view
pics = []; outDir = os.path.dirname(args.out[0])
for view in ['m','l','v','d']:
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
	#Crop medial and lateral
	if img <= 1:
		coords = (0,0,intRound(width*0.812),height)
	#Crop ventral
	elif img == 2:
		coords = (intRound(width*0.188),0,intRound(width*0.656),height)
	#Crop dorsal. Leave in scalar bar
	else:
		coords = (intRound(width*0.188),0,width,height)
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
