
%iConvmatTocsv.m
%
%
% author: steeve laquitaine
%purpose: load raw data, extract data saved in structure and save them as
%         independent mat variables, convert each mat variable to .csv files
%
%  input: 
%
%       rawpath  : raw data path (on local computer)
%       irootpath: "projInference" local path (where cloned)
%         subject: whom data it is (e.g., "sub01")
%          
%  usage:
%      
%       subject = 'sub03';
%       rawpath = '~/data/dataPsychophy/proj01_priorStrength/data'; 
%       irootpath = '~/proj/steeve/';
%       iloadRawSaveInIproj(rawpath,irootpath,subject)% 
%       iConvmatTocsv(irootpath,subject)
%
%dependencies :
%
%   the code requires "slDataMunging" library
%   > git clone https://github.com/steevelaquitaine/sldataMunging.git

function iConvmatTocsv(myrootpath,subject)

%load raw data
cd([myrootpath 'projInference/data/eye/' subject '/']) 
load(['dataEy' subject(end-1:end) '.mat']);

%sort raw data from structure to variable in 
%separate.mat files
directions = dataEy.d; save('directions','directions')
coherences = dataEy.coh; save('coherences','coherences')
priorstd = dataEy.pstd; save('priorstd','priorstd')
eyexy = [dataEy.eyexPos dataEy.eyeyPos]; save('eyexy','eyexy')

%save each separate variable in a separate .csv files
slconvMatToCsv(directions,'directions')
slconvMatToCsv(coherences,'coherences')
slconvMatToCsv(priorstd,'priorstd')
slconvMatToCsv(eyexy,'eyexy')


