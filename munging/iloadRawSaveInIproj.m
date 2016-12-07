
%iloadRawSaveInIproj.m
%
%
% author : steeve laquitaine
%purpose : load raw data from local computer and save in "projInference"
%          project
%          raw data are an .edf and .mat files containing eye positions 
%          and tasks informations
%
%
%  usage : 
%
%     subject = 'sub01';
%     rawpath = '~/data/dataPsychophy/proj01_priorStrength/data'; 
%     irootpath = '~/proj/steeve/';
%     iloadRawSaveInIproj(rawpath,irootpath,subject)
%
%
%  input :
%
%       rawpath  : raw data path (on local computer)
%       irootpath: "projInference" local path (where cloned)
%         subject: whom data it is (e.g., "sub01")

function iloadRawSaveInIproj(rawpath,irootpath,subject)

%load sorted data and save in "projInference"
dataEy = SLanalysesEyeMvt({subject},'dataPath',rawpath); 
mkdir([irootpath 'projInference/data/eye/' subject '/'])
save([irootpath 'projInference/data/eye/' subject '/dataEy' subject(end-1:end) '.mat'],'dataEy')
