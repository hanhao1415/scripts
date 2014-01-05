#!/usr/bin/python
"""

Save paired grid views of FreeSurfer surface.
Image pasting code modified from http://29a.ch/2009/5/14/concatenating-images-using-python

Note: Scale bar code still needs some work.
"""

#Parse some arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Combine Pysurfer Views')
#Positional arguments
arg_parse.add_argument('lh',help='Left hemisphere statistical map in FreeSurfer surface space.',nargs=1)
arg_parse.add_argument('rh',help='Right hemisphere statistical map in FreeSurfer surface space.',nargs=1)
arg_parse.add_argument('out',help='Name of output image.',nargs=1)
#Optional arguments
arg_parse.add_argument('-surf',help='Freesurfer surface to display map on. Default is pial.',nargs=1,default=['pial'])
arg_parse.add_argument('-minmax',help='Minumum and maximum for map. Default is 2 and 5.',nargs=2,default=[2,5],type=float)
arg_parse.add_argument('-sign',help='Sign for overlay. By default it is determined from data',nargs=1,default=[0])
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


	if ( sign == 'abs' or sign == 'pos' ):
		overlay.pos_bar.scalar_bar_representation.position2 = (0.15,0.45)
		overlay.pos_bar.scalar_bar_representation.position = (0.82,0.55)
		overlay.pos_bar.scalar_bar_representation.orientation = 1
		overlay.pos_bar.lut.table = autumnData * 255
	if ( sign == 'abs' or sign == 'neg' ):
		overlay.neg_bar.scalar_bar_representation.position2 = (0.15,0.45)
		overlay.neg_bar.scalar_bar_representation.position = (0.82,0.075)
		overlay.neg_bar.scalar_bar_representation.orientation = 1
		overlay.neg_bar.lut.table = invertData * 255

pics = []; outDir = os.path.dirname(args.out[0])
for hemi in ['lh','rh']:

	#Get brain image up
	brain = surf.Brain('fsaverage',hemi,args.surf[0],config_opts={'size':'800','background':'white'})

	#Read in data
	if hemi == 'lh':
		data = surf.io.read_scalar_data(args.lh[0])
	else:
		data = surf.io.read_scalar_data(args.rh[0])
	
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
		pic = os.path.join(outDir,'temp_%s_%s.jpg'%(hemi,view))
		brain.save_image(pic)
		pics.append(pic)

#Reorder lists
pics = [pics[0],pics[2],pics[1],pics[3]]

#Read in each image
images = map(Image.open, pics)

#Crop each image
cropped = []
for img in range(len(images)):
	coords = (0,160,600,565)
	crop = images[img].crop(coords)
	cropped.append(crop)

#Crop out the scale bar
timeImage = images[0].crop((650,25,775,375))

#Make a new output image with a white background
output = Image.new("RGB",(1375,825))
output.paste((255,255,255))

#Paste everything to new image
output.paste(cropped[0],(0,5))
output.paste(cropped[1],(600,0))
output.paste(cropped[2],(0,405))
output.paste(cropped[3],(600,405))
output.paste(timeImage,(1215,300))

#Save the grid figure
output.save(args.out[0])

#Remove the temp files
for pic in pics:
	os.remove(pic)