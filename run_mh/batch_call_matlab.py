#!/bin/bash
# Author: Brian Munsky
# Date: 09/24/16
# Purpose: Batch Matlab Job Submission

import subprocess
import time
for i in [0,1,2,3,4,5,6,7,8,9,10]:  #population size
	cmd = ' '.join( ['qsub','-q munsky.q@node* -cwd -o /dev/null -e Errs/','-v','','num=%i' % (i),'Matlab_Wrapper.sh'] )
	print cmd
	subprocess.call( cmd, shell=True )
	time.sleep(0.1) # delay for 0.1 seconds.	
