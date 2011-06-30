#!/usr/bin/python
import csv,sys

#import motion paramaters file and skip the header row
motion_file = csv.reader(open(sys.argv[1],'r'),delimiter=' ',quotechar="'",quoting=csv.QUOTE_NONNUMERIC)
motion_file.next()

#put each column into a list
frame,disx,disy,disz,rotx,roty,rotz = zip(*motion_file)
#loop through every frame, checking for translation and rotation outliers
dis_outliers = 0; rot_outliers = 0
for x_dis,y_dis,z_dis,x_rot,y_rot,z_rot in zip(disx,disy,disz,rotz,roty,rotz):
	if abs(float(x_dis))>=1 or abs(float(y_dis))>=1 or abs(float(z_dis))>=1:
		dis_outliers += 1
	if abs(float(x_rot))>=.5 or abs(float(y_rot))>=.5 or abs(float(z_rot))>=.5:
		rot_outliers += 1

print dis_outliers,rot_outliers
