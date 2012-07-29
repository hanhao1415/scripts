#!/data/nil-bluearc/benzinger2/Tyler/epd-7.2-2-rh5-x86_64/bin/python

###########
###Usage###
###########
import argparse
arg_parse = argparse.ArgumentParser(description='Eta Generator')
#Positional Arguments
arg_parse.add_argument('params',help='Param file with location of all images',nargs=1)
arg_parse.add_argument('weights',help='Weight file with weights',nargs=1)
args = arg_parse.parse_args()

#Import necessary libraries
import re, numpy as np, nibabel as nib, sys, subprocess as sp, os, shutil

#Load in params file
try:
	execfile(args.params[0])
except (IOError):
	print 'Cannot find parameter file. Exiting...'
	sys.exit()

#Read in weights file
try:
	with open(args.weights[0],'r') as file:
    		weights = file.read().splitlines()
except IOError:
	print 'Cannot find weight file. Exiting...'
	sys.exit()

#Strip away all characters that are not -,., or a number
for weight in range(len(weights)):
	weights[weight] = re.sub("[^0-9.-]", "",weights[weight])

#Check to see there are enough weights. Will probably fail with blanks.
numWeights = int(weights[1])
if (len(weights) < 2+numWeights):
	print 'Number of weights is less than indicated in weight file. Exiting...'
	sys.exit()

#Remove uncessary rows
weights = [weights[0]]+weights[2:2+numWeights]

#Load in basis library and mean image
basis = nib.load(basis); basisData = basis.get_data()
mean = nib.load(mean); meanData= mean.get_data()

#Make field map image using the current set of weights
meanShape = meanData.shape
scaledData = np.sum(basisData*np.array(weights[1::],dtype=np.float32),axis=3) + \
	meanData.reshape(meanShape[0],meanShape[1],meanShape[2])
scaled = nib.Nifti1Image(scaledData,mean.get_affine())
scaled.to_filename('scaled_phase.nii')

#Create log file
log = open('basis_opt.log','w')

#Set output type to nii to appease nifti_4dfp
os.putenv('FSLOUTPUTTYPE','NIFTI')

#Shift median to 0
stat = sp.Popen('fslstats scaled_phase -k %s -P 50'%(mask),shell=True,stdout=sp.PIPE,stderr=log)
median = float(stat.communicate()[0].rstrip())
sp.call('fslmaths scaled_phase -sub %s -mas %s scaled_phase'%(median,mask),\
	shell=True,stdout=log,stderr=log)

#Run fugue with the new field map
sp.call('fugue --loadfmap=scaled_phase --dwell=%s -u epi_unwarped -i %s --mask=%s --unwarpdir=%s'\
	%(float(weights[0])/1000,epi,mask,dir),shell=True,stdout=log,stderr=log)

#Convert undistorted epi to 4dfp
sp.call('nifti_4dfp -4 epi_unwarped epi_unwarped -N',shell=True,stdout=log,stderr=log)

#Copy over the initial t4
if os.path.isfile('epi_unwarped_to_t2_t4'):
	os.remove('epi_unwarped_to_t2_t4')
try:
	shutil.copy(t4,'epi_unwarped_to_t2_t4')
except IOError:
	print 'Cannot copy t4 file. Perhaps you do not have write permissions?'
	sys.exit()
				
#Register the unwarped epi to the highres
sp.call('imgreg_4dfp %s none epi_unwarped none epi_unwarped_to_t2_t4 8195 > \
	epi_unwarped_to_t2.log'%(t2),shell=True,stdout=log,stderr=log)

#Extract eta from imgreg log
awk = sp.Popen('awk \'/eta,q/{print $2}\' epi_unwarped_to_t2.log | tail -1',shell=True,\
	stdout=sp.PIPE,stderr=log)
eta = float(awk.communicate()[0].rstrip())

#Recreate weight file and append eta to the end of it
file = open(args.weights[0],'w')
file.write('Echo Spacing: %.5f'%(float(weights[0])))
file.write('\nWeights: %i'%(numWeights))
for weight in range(numWeights):
	file.write('\n %15.5f'%(float(weights[weight+1])))
file.write('\nEta: %5f'%(eta))
file.close()
