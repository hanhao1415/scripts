#!/usr/bin/python 

#Parse some arguments
import argparse
arg_parse = argparse.ArgumentParser(description='Combine Medial and Lateral Movies Snapshots')
#Positional arguments
arg_parse.add_argument('med',help='Directory containing medial snapshots',nargs=1)
arg_parse.add_argument('lat',help='Directory containing lateral snapshots',nargs=1)
arg_parse.add_argument('out',help='Root for outputed images',nargs=1)
arg_parse.add_argument('-pib',action='store_const',const=1,help='Crop using PiB presets')
args = arg_parse.parse_args()

#Load the necessary libraries
import Image, os

#Get a list of all the images in the directory
medImages = os.listdir(args.med[0])
latImages = os.listdir(args.lat[0])

#Remove the DS_STORE files 
try:
	medImages.remove('.DS_Store')
	latImages.remove('.DS_Store')
except ValueError:
	print '.DS_Store files are not present in both lists. This should be ok...'

for med,lat in zip(medImages,latImages):
	
	outImage = "%s_%s.jpeg"%(args.out[0],os.path.splitext(med)[0])
	print 'Making image %s'%(outImage)

	if (med != lat):
		print 'Error: Medial and lateral images do not have the same file name. Exiting...'
		sys.exit()
	
	#Make a new combined image
	combined = Image.new("RGB",(1366,768))

	#Load in two images
	images = map(Image.open,[os.path.join(args.med[0],med),os.path.join(args.lat[0],lat)])
	
	#Crop and paste images
	medImage = images[0].crop((40,160,610,575))
	latImage = images[1].crop((40,160,610,575))
	if ( args.pib == 1 ):
		scaleImage = images[0].crop((650,50,800,550))
		timeImage = images[0].crop((90,600,715,675))
	else:
		scaleImage = images[0].crop((650,0,800,750))
		timeImage = images[0].crop((40,600,650,675))

	#Paste croped images to combined image
	combined.paste(medImage,(30,145))
	combined.paste(latImage,(580,145))
	combined.paste(scaleImage,(1180,5))
	combined.paste(timeImage,(315,615))
	
	#Write out new image
	combined.save(outImage)