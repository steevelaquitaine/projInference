

%SLgetActivemFile.m
%
%author: steeve laquitaine
%  date: 140801
% purpose: get active m file

function filename = SLgetActivemFile

filename = matlab.desktop.editor.getActiveFilename;
