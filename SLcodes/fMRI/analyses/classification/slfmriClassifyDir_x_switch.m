
%slfmriClassifyDir_x_switch.m
%
%
% author: steeve laquitaine
%purpose: classify motion directions from teh brain activity for data 
%         sorted by "switching" variable to determine whether sensory 
%         evidence is represented or not when subject switch to prior.
%           
%
%  usage : 
%
%       e.g., 
%
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%           [dataClassifsw1,dataClassifsw2] = slfmriClassifyDir_x_switch(o)
%
%
%       e.g., 
%
%           %init data info, set analysis parameters and stack instances
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%           params = {'accuracyAtTime',[7 14],'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1'};           
%           myClass = {{'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06'},....
%                      {'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06'}};
% 
%           o1 = slfmriGetClassifAnalysisParams(o,myClass{1},params{:});
%           [o1,~,ssw1] = slStackSessionInstances(o1,params{:});
%           
%           o2 = slfmriGetClassifAnalysisParams(o,myClass{2},params{:});
%           [o2,~,ssw2] = slStackSessionInstances(o2,params{:});
%           


function [dataClassifsw1,dataClassifsw2] = slfmriClassifyDir_x_switch(o)

%weak motion (6 % coherence)

%--------- switch to prior -------

%Dataset was balanced between classes
myClass = {'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06'};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

%load results
cd([o.classifResultPath '/' cell2mat(myClass) '/accAtTimeleaveOneOutfisherbalancByRemovI/' o.myROIname{1}])
dataClassifsw1 = load(['ClassifStckSess' cell2mat(myClass) '_' o.myROIname{1} '.mat']);



%--------- switch to motion direction -------

%Dataset was balanced between classes
myClass = {'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06'};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

%load results
cd([o.classifResultPath '/' cell2mat(myClass) '/accAtTimeleaveOneOutfisherbalancByRemovI/' o.myROIname{1}])
dataClassifsw2 = load(['ClassifStckSess' cell2mat(myClass) '_' o.myROIname{1} '.mat']);

