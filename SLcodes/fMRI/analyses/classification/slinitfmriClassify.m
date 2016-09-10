

%slinitfmriClassify.m
%
%    Author: steeve laquitaine
%      date: 151208
%   purpose: initialize classification for fMRI data (sets session and params)
%
%inputs
%
%       'stckSessions'
%
%in progress ...


function slinitfmriClassify(varargin)

%get args
getArgs(varargin,{'gr','sc','roi','tnum','pnum','snum','rtp','ses'});

keyboard

%init stack sessions
if slIsInput(varargin,'stckSessions')
    o.myGroup = 'Concatenation';
    o.rootpath = '/Volumes/Transcend/data/sltaskdotdirfmri05/';
    o.sessPath{1} = [o.rootpath '/s02520150814'];
    o.sessPath{2} = [o.rootpath '/s02520150814'];
    o.sessPath{3} = [o.rootpath '/s02520150814'];
    o.stckPath = o.rootpath;
    o.myScan = 1;
    o.myROIname = {'lV1'};%,'V3A','MT','IPS'};
    o.taskNum = 2;
    o.phaseNum = 1;
    o.segmentNum = 2;  
end
