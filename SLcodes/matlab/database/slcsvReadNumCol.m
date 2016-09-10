

%slcsvReadNumCol.m
%
%
%author : steeve laquitaine
%  date : 160714
%purpose: count the number of columns in a csv file
%
%
% usage:
%
%       nCols =  slcsvReadNumCol('modelfitData.csv')

function nCols = slcsvReadNumCol(filename)

%First we need to get the header-line
fid1 = fopen(filename, 'r');
Header = fgetl(fid1);
fclose(fid1);

%Convert Header to cell array
Header = regexp(Header, '([^,]*)', 'tokens');
Header = cat(2, Header{:});
nCols = length(Header);