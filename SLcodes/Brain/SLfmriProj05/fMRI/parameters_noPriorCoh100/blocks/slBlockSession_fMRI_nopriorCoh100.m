
%slBlockSession_fMRI_nopriorCoh100.m
%
% author: steeve laquitaine	
%   date: 160423
%purpose: divide prior with mean 135 deg fMRI session into N blocks of
%         15 trials (and a short last block)
%  usage: 
%	
%			slBlockSession_fMRI_nopriorCoh100


%set file parameter name
param_nm = 'fMRIparamsNoPriorCoh100';



%% load parameters
load([param_nm '.mat'])

taskAll = task;
clear task

%randomize sequences
myRand = randperm(taskAll.parameter.dir.trialnum);
taskAll.parameter.dir.series = taskAll.parameter.dir.series(myRand);
taskAll.parameter.dir.coh = taskAll.parameter.dir.coh(myRand);

%make N blocks of 15 trials each
nTrialPerBlock = 15;
nBlocks = round(taskAll.parameter.dir.trialnum/nTrialPerBlock);
blocki  = repmat([1:nBlocks]',1,15)';
lastBlock = [];

if length(blocki(:)) ~= taskAll.parameter.dir.trialnum    
    nTrialLastBlock = taskAll.parameter.dir.trialnum - length(blocki(:));
    lastBlock = repmat(nBlocks+1,1,nTrialLastBlock);
    nBlocks = nBlocks + 1;
end

blocks = [blocki(:); lastBlock'];
blocksi = unique(blocks);

%save to 'block' folder
mypath = which([param_nm '.mat'])
mypath = SLgetFileParentDir(SLgetFileParentDir(mypath));
cd([mypath '/blocks'])

%create blocksi
for i = 1 : nBlocks    
    clear task
    thisBlock = blocks == blocksi(i);    
    task.parameter.dir.series = taskAll.parameter.dir.series(thisBlock);
    task.parameter.dir.coh = taskAll.parameter.dir.coh(thisBlock);
	task.parameter.dir.mean = taskAll.parameter.dir.mean;
	task.parameter.dir.std = taskAll.parameter.dir.std;
	task.parameter.dir.sample = unique(taskAll.parameter.dir.series);
	task.parameter.dir.block = i;

	%name
	if i >= 10
		fname = [param_nm '_block' num2str(i) '.mat'];
	else
	 	fname = [param_nm '_block0' num2str(i)];
	end
	save(fname,'task','taskAll');    
end
fprintf('%s \n','done')