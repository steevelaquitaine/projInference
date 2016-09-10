

%slfmriLoadInstancedb.m
%
%
% author: steeve laquitaine
%   date: 160630
% pupose: load instance database
%
%usage: 
%
%       slfmriInitAnalysisTaskDotDirfMRI05
%       [o,behbySess] = slmakeSwitchingVar(o);  
%       [~,o] = slfMRImakeDatabase(o,varargin);
%       d = slfmriLoadInstancedb('loadInstancesdb')

function d = slfmriLoadInstancedb(varargin)

if any(strcmp(varargin,'loadInstancesdb'))
    d = importdata('d.mat');
elseif any(strcmp(varargin,'useinputdb'))
    fprintf('%s \n','(slfmriGetVoxPopActivity2) Using input database.')
    d = varargin{find(strcmp(varargin,'useinputdb'))+1};
end