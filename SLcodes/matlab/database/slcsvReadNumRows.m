

%slcsvReadNumCol.m
%
%
%author : steeve laquitaine
%  date : 160714
%purpose: count the number of rows in a csv file
%
%
% usage:
%
%       nRows =  slcsvReadNumRows('modelfitData.csv')

function nRows = slcsvReadNumRows(filename)

fid = fopen(filename);
allText = textscan(fid,'%s','delimiter','\n');
nRows = length(allText{1});
fclose(fid)

