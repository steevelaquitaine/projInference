


%slfmriCreateOrLoadfMRIdb.m
%
%
% author : steeve laquitaine
%purpose : create fMRI database "d" from fMRI tSeries and .mat files
%          or load existing database structure "d" from workspace
%
%   usage: 
% 
%           slfmriInitAnalysisTaskDotDirfMRI05;                      
%           load d.mat           
%           [d,o]= slfmriCreateOrLoadfMRIdb(o,[varargin,'db',d])


function [d,o]= slfmriCreateOrLoadfMRIdb(o,varargin)

%get the database over sessions
if ~any(strcmp(varargin{1},'db'))
    [d,o] = slfmriGetDBoverSessions(o,varargin);
else
    d = varargin{1}{find(strcmp(varargin{1},'db')==1)+1};
end
