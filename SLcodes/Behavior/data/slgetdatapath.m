
%slgetdatapath.m
%
%
% author: steeve laquitaine
%purpose: get the data path for the motion direction prior project
%
%
% usage :
%
%           datapath = slgetdatapath('lab')
%           datapath = slgetdatapath('home')


function  datapath = slgetdatapath(where)

if strcmp(where,'lab')
    datapath = '~/data/dataPsychophy/proj01_priorStrength/data/';
end
if strcmp(where,'home')
    datapath = '/Dropbox/myDropbox/Project_withJustin/data/dataPsychophy/proj01_priorStrength/data/';
end
