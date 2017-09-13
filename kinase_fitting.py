#!/bin/sh
# Author: Brian Munsky/Zach Fox
# Date: 10/14/14
# Purpose: Batch Python Job Submission

import subprocess
import time

minjob = 1

maxjob = 50 
filenames ='submit_full_fit'
path = 'out/'
#print filenames
for taskID in range(minjob,maxjob):
        print taskID
        cmd = ' '.join( ['qsub','-q','munsky.q@node*','-o',path+filenames+'_%d.dat -N' % taskID,filenames+'_%d' % taskID,'-e',path+filenames+'_%d.error' %taskID,' -v taskID=%d ' % taskID,filenames+'.sh','',str(taskID)] )
        #The command is "qsub SEND_VARIABLE=X script.sh"
        
	subprocess.call( cmd, shell=True )
	time.sleep(2) # delay for 0.1 seconds for random number generator seed LOL
