
%slfmriStackAndSaveSessionInstances.m
%
% author: steeve laquitaine
%purpose: stack and save session instances in a structured directory
%
%  usage :
%
%
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%           params = {'accuracyAtTime',[7 14],'loadSavedROI','CalcInstances',...
%                     'leaveOneOut','fisher','balancByRemovI=1',...
%                     'numInstances',8,'loadInstancesByCond'};
%           myConds = {{'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06'},....
%                      {'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06'}};
%           [cbyc,sbyc,o] = slfmriStackAndSaveSessionInstances(o,myConds,params)


function [obyc,sbyc,o] = slfmriStackAndSaveSessionInstances(o,myConds,params)

%get stacked instances over sessions for each condition
for myc = 1 : length(myConds)
    o = slfmriGetClassifAnalysisParams(o,myConds{myc},params{:});
    [obyc{myc},~,sbyc{myc}] = slStackSessionInstances(o,params{:});
end
%save instances by conditions
o.analysispath = ([o.rootpath '/slStckAnalyses/Concatenation/classif/slfmriClassifyByConditions/' o.myROIname{1} '/' cell2mat([myConds{:};])  '/'  ]);
mkdir(o.analysispath); cd(o.analysispath)
save('dataInstancesByCond','obyc','sbyc','o','params','myConds')
