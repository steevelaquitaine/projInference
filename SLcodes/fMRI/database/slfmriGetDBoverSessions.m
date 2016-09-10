%slfmriGetDBoverSessions.m
%
%
% author: steeve laquitaine
%   date: 160729
%purpose: get database of :
%         - trial instances of voxel response to the motion stimulus
%         - associated variables (motion direction, coherence, subject switching to prior or evidence)        
%
%usage : 
%           slfmriInitAnalysisTaskDotDirfMRI05
%           [d,o,behbySess] = slfmriGetDBoverSessions(o,varargin)

function [d,o,behbySess] = slfmriGetDBoverSessions(o,varargin)

%% calculate/add "switch" variable to sessions (once)
[o,behbySess] = slmakeSwitchingVar(o);  

%% make database
[d,o] = slfMRImakeDatabase(o,varargin{:});
