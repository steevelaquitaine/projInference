


%slPrintTxt.m
%
%
% author: steeve laquitaine
%purpose: print variable x number value in a text file with name filename
%
%
%   usage: 
%
%       slPrintTxt('mytext.txt',rand(10))


function slPrintTxt(filename,x)

fileID = fopen(filename,'w');
fprintf(fileID,'%.2f ',x);
edit(filename)