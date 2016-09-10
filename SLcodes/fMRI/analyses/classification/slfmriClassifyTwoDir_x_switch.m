
%slfmriClassifyTwoDir_x_switch.m
%
%author: steeve laquitaine
%purpose: classify two motion directions 15 and 85 deg from brain activity 
%         for data sorted by "switching" variable to determine whether 
%         sensory evidence is represented or not when subject switch to prior.
%           
%
%  usage : 
%
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%           [dataClassifsw1,dataClassifsw2] = slfmriClassifyTwoDir_x_switch(o)

function [dataClassifsw1,dataClassifsw2] = slfmriClassifyTwoDir_x_switch(o)

%weak motion (6 % coherence)

%--------- switch to prior -------

%Dataset was balanced between classes
myClass = {'myRandomDir=15_x_mySwitch=1_x_myRandomCoh=0.06','myRandomDir=85_x_mySwitch=1_x_myRandomCoh=0.06'};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

% slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
%     'CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

%load results
cd([o.classifResultPath '/' cell2mat(myClass) '/accAtTimeleaveOneOutfisherbalancByRemovI/' o.myROIname{1}])
dataClassifsw1 = load(['ClassifStckSess' cell2mat(myClass) '_' o.myROIname{1} '.mat']);



%--------- switch to motion direction -------

%Dataset was balanced between classes
myClass = {'myRandomDir=15_x_mySwitch=2_x_myRandomCoh=0.06','myRandomDir=85_x_mySwitch=2_x_myRandomCoh=0.06'};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

% slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
%     'CalcInstances','leaveOneOut','fisher','balancByRemovI=1');


%load results
cd([o.classifResultPath '/' cell2mat(myClass) '/accAtTimeleaveOneOutfisherbalancByRemovI/' o.myROIname{1}])
dataClassifsw2 = load(['ClassifStckSess' cell2mat(myClass) '_' o.myROIname{1} '.mat']);

