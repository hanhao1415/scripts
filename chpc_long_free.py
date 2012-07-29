#!/usr/bin/python

import sys,os,csv,argparse,textwrap,subprocess as sp

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Setup Longitudinal FreeSurfer Batch Scripts')
arg_parse.add_argument('dir',help='Directory containing Freesurfer sessions',nargs=1)
arg_parse.add_argument('csv',help='CSV file containing subject ID and session IDs',nargs=1)
arg_parse.add_argument('-wtime',help='Walltime in hours for each job. Default is 48',nargs=1,type=int,default=[48])
arg_parse.add_argument('-vmem',help='Amount of virtual memory in GB for job. Default is 4',nargs=1,type=int,default=[4])
arg_parse.add_argument('-aff',help='Use affine registration to the base. Default is rigid.',nargs='?',const='-base-affine',default='-base')
args = arg_parse.parse_args()

#Read CSV File
try:
	csv_reader = csv.reader(open(args.csv[0],'r'))
except (IOError):
	print 'Error: Cannot find CSV file at %s'%(args.csv[0])
	sys.exit()

#Check if SUBJECTS_DIR exists
subj_dir = os.path.abspath(args.dir[0])
if os.path.exists(subj_dir) is False:
	print 'Error: Cannot find directory at %s'%(args.dir[0])
	sys.exit()
	
#Make directory for template batch files
long_dir = os.path.join(subj_dir,'long'); batch_dir = os.path.join(long_dir,'batch')
log_dir = os.path.join(long_dir,'logs')
try:
	os.mkdir(long_dir); os.mkdir(batch_dir); os.mkdir(log_dir)
except (OSError):
	print 'Long directory already exists in %s. This is ok, but present files may be overwritten.'%(subj_dir)

#Loop through each line/subject
for item in csv_reader:
	
	#Skip blank lines
	if len(item) == 0: continue
	
	#Extract subject
	subj = item[0]; item.remove(subj);
			
	#Create batch file for each time point
	tp_string = ''; ok = 0; tp_paths = []
	for tp in item:
		if ok== 1: continue
	
		#Check for session's existence
		if os.path.exists(os.path.join(subj_dir,tp)) is False:
			print 'Error: Cannot find session %s with %s. Skipping subject %s.'%(tp,subj_dir,subj)
			ok = 1; continue 
		
		#Setup the time point string for later batch processing
		tp_string = '%s -tp %s'%(tp_string,tp)
	
		#Create the batch file
		tp_batch = textwrap.dedent('''\
			#Longitudinal Batch Script for Timepoint: %s of Base: %s
			
			#Setup PBS options
			#PBS -N %s -l nodes=1:ppn=1:idataplex,walltime=%i:00:00,vmem=%igb -o %s/%s_stdout.txt -e %s/%s_stderr.txt
			
			#Setup FreeSurfer options
			export SUBJECTS_DIR=%s
			export FREESURFER_HOME=/export/freesurfer-5.1
			export PATH=$FREESURFER_HOME/bin/:$PATH
			export PATH=${FREESURFER_HOME}/mni/bin:$PATH
			export PERL5LIB=${FREESURFER_HOME}/mni/lib/perl5/5.8.5/

			#Run the command
			recon-all -long %s %s -qcache -all
			'''%(tp,subj,tp,args.wtime[0],args.vmem[0],log_dir,tp,log_dir,tp,subj_dir,tp,subj))
		
		#Write out batch file
		tp_path = os.path.join(batch_dir,'%s_batch.txt'%(tp)); tp_paths.append(tp_path)
		file = open(tp_path,'w+'); file.write(tp_batch); file.close()
	
	if ok == 1: continue
	
	#Create the batch file for the template
	temp_batch = textwrap.dedent('''\
			#Longitudinal Batch Script for Base: %s
			
			#Setup PBS options
			#PBS -N %s -l nodes=1:ppn=1:idataplex,walltime=%i:00:00,vmem=%igb -o %s/%s_stdout.txt -e %s/%s_stderr.txt

			#Setup FreeSurfer options
			export SUBJECTS_DIR=%s
			export FREESURFER_HOME=/export/freesurfer-5.1
			export PATH=$FREESURFER_HOME/bin/:$PATH
			export PATH=${FREESURFER_HOME}/mni/bin:$PATH
			export PERL5LIB=${FREESURFER_HOME}/mni/lib/perl5/5.8.5/

			#Run the command
			recon-all %s %s %s -qcache -all
			'''%(subj,subj,args.wtime[0],args.vmem[0],log_dir,subj,log_dir,subj,subj_dir,args.aff,subj,tp_string))
	base_file = os.path.join(batch_dir,'%s_batch.txt'%(subj))
	file = open(base_file,'w+'); file.write(temp_batch); file.close()
	
	#Submit the template and save the job id
	print 'Launching jobs for subject: %s'%(subj)
	qsub = sp.Popen('qsub %s'%(base_file),shell=True,stdout=sp.PIPE)
	jid = qsub.communicate()[0].rstrip()
	
	#Submit the timepoints. Hold for the template job ID.
	for path in tp_paths:
		qsub = sp.call('qsub -W depend=afterok:%s %s'%(jid,path),shell=True,stdout=sp.PIPE)
print 'All jobs launched....'






