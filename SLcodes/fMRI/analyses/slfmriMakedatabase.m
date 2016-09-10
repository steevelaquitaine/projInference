
%slfmriMakedatabase.m
%
% author: steeve laquitaine
%   date: 160108
%purpose: create database matrix of factors and data
%
%%usage: 
%
%       slfmriMakedatabase(varname,taskNum,phaseNum,segmentNum,sessPath)
%
%
%
%       sessPath = '/Volumes/Transcend/data/sltaskdotdirfmri05/';
% 
%       slfmriMakedatabase({'myRandomCoh','myRandomDir'},'taskNum=2',...
%               'phaseNum=1','segmentNum=2',sessPath)
%
%option input:
%       'stacked': data are stacked over sessions

%
%in progress

function slfmriMakedatabase(varargin)

%getArgs
o.nPerm = [];
getArgs(varargin{'varname','accuracyAtTime','taskNum','phaseNum',...
    'segmentNum','sessPath'});

%get classes
o.varname = varname;
o.taskNum = varname;
o.phaseNum = varname;
o.segmentNum = varname;
o.sessPath = sessPath;

%load session
cd(o.sessPath)
h = mrLoadRet([]);
fprintf('\n %s \n \n', ['(slCalcStackedInstances) Session:' o.sessPath{si} ' ------'])

%get stimvols by class
for ivar = 1 : length(varname)
    
    [d{ivar},o.vars{ivar}] = getStimvol(getMLRView,o.varname,'taskNum',o.taskNum,...
        'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
    
end

dbstack
keyboard