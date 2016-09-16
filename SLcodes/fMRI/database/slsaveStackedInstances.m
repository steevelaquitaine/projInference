
%slsaveOrLoadStackedInstances.m
%
%
%
% author : steeve laquitaine
%purpose : save instances from stacked sessions



function o = slsaveOrLoadStackedInstances(o,varargin)

%path
o.savedInstancesPath = [o.rootpath 'slStckAnalyses/instanceMatrix/' o.sid '/' o.myROIname '/'];

%move to path
cd(o.savedInstancesPath)

%save
save('instanceMatrix','d')

