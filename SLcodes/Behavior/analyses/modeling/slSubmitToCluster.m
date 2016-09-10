
%slSubmitToCluster.m
%
%
% author: steeve laquitaine
%   date: 151022
% status: Incomplete
%purpose: run matlab command in myfun on Stanford Sherlock cluster
%
% usage:
%
%       slSubmitToCluster('slFitWJMlgTprior','sub12','[1 20 20 20 20 0.1 60 0.4]','167:00:00','long','1','16')
%       slSubmitToCluster('testSherlock','JustTesting','[1 20 20 20 20 0.1 60 0.4]','167:00:00','long','1','16')
%
%note:
%   
%     To do : Does not seem to work yet because ' and " characters are apparently 
%     not in the correct format.....
%     
%
%make sure your code library is updated on the cluster


function [fname,cmdout3] = slSubmitToCluster(myfun,subject,initP,mytime,qos,nodes,ntaskspernode)

%batch file
fname = [subject '.matlab'];

%remove existing submission file
cd ~/Desktop
[status0,cmdout0] = system(sprintf('%s','rm ',fname));

%write batch file
system(sprintf('%s %s','echo "#!/bin/bash" >>',fname))
system(sprintf('%s%s','echo "#SBATCH --job-name=',subject,'">>',fname))
system(sprintf('%s%s','echo "#SBATCH --output=',subject,'.out">>',fname))
system(sprintf('%s%s','echo "#SBATCH --error=',subject,'.err">>',fname))
system(sprintf('%s%s','echo "#SBATCH --time=',mytime,'">>',fname))
system(sprintf('%s%s','echo "#SBATCH --qos=',qos,'">>',fname))
system(sprintf('%s%s','echo "#SBATCH --nodes=',nodes,'">>',fname))
system(sprintf('%s%s','echo "#SBATCH --ntasks-per-node=',ntaskspernode,'">>',fname))
system(sprintf('%s%s','echo "#SBATCH --export=All">>',fname))
system(sprintf('%s%s','echo "#SBATCH --mem-per-cpu 4000">>',fname))
system(sprintf('%s%s','echo "#SBATCH --mem 64000">>',fname))
system(sprintf('%s%s','echo "#SBATCH --exclusive">>',fname))
system(sprintf('%s%s','echo "#SBATCH --mail-type=ALL">>',fname))

system(sprintf('%s%s','echo "##my commands">>',fname))
system(sprintf('%s%s','echo "##load matlab">>',fname))
system(sprintf('%s%s','echo "ml load matlab">>',fname))

system(sprintf('%s%s','echo "##add code library path">>',fname))
strsub = ['''' subject ''''];
system(sprintf('%s%s%s','echo "matlab -nojvm -r \"','addpath(genpath(',['''' '~/proj' ''''],')',...
    ');',myfun,'({',strsub,'},',initP,',',['''' 'experiment' ''''],...
    ',',['''' 'vonMisesPrior' ''''],',',['''' 'sherlock' '''',')\" output',subject,'.log">>',fname]))

%remove existing file in cluster
[status1,cmdout1] = system(sprintf('%s','ssh steeve@sherlock rm ',fname))

%send to cluster
command = sprintf('%s%s %s','rsync -a --stats ~/Desktop/',fname,'steeve@sherlock.stanford.edu:');
fprintf('%s \n',sprintf('(slSubmitfitToSherlock) Doing remote command: %s',command));
[status2,cmdout2] = system(command,'-echo')


%submit job
command = sprintf('%s','ssh steeve@sherlock sbatch ',fname);
fprintf('%s \n',sprintf('(slSubmitfitToSherlock) Doing remote command: %s',command));
[status3,cmdout3] = system(command,'-echo')



