#!/bin/bash
#
#/usr/local/MATLAB/R2015b/bin/matlab -nosplash -nodisplay -r "i=6;j=3;Batch_Matlab_Example(i,j,'Results/Output');exit;" > Logs/Log_6_3.txt
/usr/local/MATLAB/R2020a/bin/matlab -nosplash -nodisplay -r "i=${i},run_mh_cluster(i);exit;" > Logs/Log_${i}.txt
