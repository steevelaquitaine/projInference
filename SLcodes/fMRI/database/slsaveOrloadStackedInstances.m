
%slsaveOrloadStackedInstances.m
%
%
%
% author : steeve laquitaine
%purpose : save calculated or load pre-existing instances from 
%          stacked sessions
%
%   usage :
%
%       o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';
%       o.sid = 's025';
%       o.myROIname = 'V1';
%       o.priormean = '225';
%       o = slsaveOrloadStackedInstances(o,'save')
%       o = slsaveOrloadStackedInstances(o,'load')

function o = slsaveOrloadStackedInstances(o,saveOrload)

%path
o.savedInstancesPath = [o.rootpath 'slStckAnalyses/instanceMatrix/' o.sid '/' o.priormean '/' o.myROIname '/'];

%move to path
cd(o.savedInstancesPath)

%save or load
if strcmp(saveOrload,'save')
    %save
    save('instanceMatrix','d')
elseif strcmp(saveOrload,'save')
    %load
    load('instanceMatrix.mat')
end

