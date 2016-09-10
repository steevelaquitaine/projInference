
%slfMRImakeDatabase.m
%
%     author: steeve laquitaine
%    purpose: create instance database with data and variables stacked over sessions
%      usage: 
%           
%           slfmriInitAnalysisTaskDotDirfMRI05
%           [d,o] = slfMRImakeDatabase(o,varargin)
%
%Description :  first initialize data info and analysis, then run function

function [d,o] = slfMRImakeDatabase(o,varargin)
tic
%extract varargin
while size(varargin)==1
    varargin = varargin{:};
end

%get database (responses and variables)
slfmriClassify2(o,o.myVar,varargin{:});

%create stacked database of instances and associated variables
%[o,sessPath] = slfmriInitGetInstancedb(o.myROIname,['ClassifStckSess' o.myVar '_' o.myROIname '.mat']);
[o,d] = slfmriStckGetInstancedb(o);
toc