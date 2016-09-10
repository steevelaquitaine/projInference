%SLtaskLocEst.m
%
% author: Steeve Laquitaine
%   date: 140506
%purpose: check if mat file exist in directory
%
% usage:
%
%       o = SLisExist('myFile')

function o = SLisExist(myFileInDir)

d = dir([myFileInDir,'.mat']);

o = isempty(d);