

%SLinitRunUniPriorLocTask.m
%
%  author: steeve laquitaine
%    date: 2015/05/01
% purpose: get a file parent directory
%   usage:
%
%			myparentdir = SLgetFileParentDir('~/proj/steeve/myfile.mat')


function myparentdir = SLgetFileParentDir(myPath)

myparentdir = fileparts(myPath);