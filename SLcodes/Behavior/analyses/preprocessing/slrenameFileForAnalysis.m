
%slrenameFileForAnalysis.m
%
%
% author: steeve laquitaine
%purpose:  rename raw stim files for routine analysis
%  usage: 
%
%
%      filepath=slrenameFileForAnalysis('140425_stim01','sub=sub01','session=08','run=37','pstd=80','pmean=225','coh=[006,012]','exp=outscanner','date=160811')



function filepath = slrenameFileForAnalysis(filename,sub,ses,run,pstd,pmean,coh,exp,dat)

%cleanup workspace
%clearex('filename','sb','ses','run','pstd','pmean','coh','date')

%load file
load(filename)

%get file info
ex = exp(find(exp=='=')+1 : end);
sb = sub(find(sub=='=')+1 : end);
ss = ses(find(ses=='=')+1 : end);
rn = run(find(run=='=')+1 : end);
ps = pstd(find(pstd=='=')+1 : end);
pm = pmean(find(pmean=='=')+1 : end);
co = coh(find(coh=='=')+1 : end); co(co==',')='_';co(co=='[')='_';co(co==']')='_';
dt = dat(find(dat=='=')+1 : end);

%...and rename
savedFile = ['steeve_fmri_data_' sb '_sess' ss '_run' rn '_Pstd0' ps ...
    '_mean' pm '_coh' co 'dir5_randInitPosSymPrior_' dt '_' ex '.mat'];

save(savedFile,'fixStimulus', 'myscreen','stimulus','task')
fprintf('%s \n','(slrenameFileForAnalysis) File was renamed for routine analysis')

filepath = which(savedFile);
