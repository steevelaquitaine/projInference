
%slcsvRead.m
%
%
% author: steeve laquitaine
%   date: 160714
%purpose: read csv table
%
%usage:
%
%       [rowheader,colheader,data] = slcsvRead('modelfitData.csv');
%
%
%Inputs : 
%
% .csv strcutured with row header, column header and data matrix :
%
%           text | text | text | ....
%           _____ ______ ______
%
%           text | num  | num  |
%           _____ ______ ______ 
%
%           text | num  | num  |
%             .
%             .
%             .

function [rowheader,colheader,data] = slcsvRead(filename)

%%%%%%%%%%%
% Headers %
%%%%%%%%%%%
%First we need to get the header-line
fid1 = fopen(filename, 'r');
colheader = fgetl(fid1);
fclose(fid1);

%Convert header to cell array
colheader = regexp(colheader, '([^,]*)', 'tokens');
colheader = cat(2,colheader{:});


%%%%%%%%%%%
% table   %
%%%%%%%%%%%
%# Read in the data
nCols =  slcsvReadNumCol(filename);
fid1 = fopen(filename, 'r');
format = ['%s' repmat('%f',1,nCols) '%s'];
data = textscan(fid1,format,'Delimiter',',','HeaderLines',1);
fclose(fid1);

%get column header
rowheader = data{1};

%convert cell column data to matrix
tmp = [];
for i = 1 : nCols - 1
    tmp = [tmp data{i+1}];
end
data = tmp;
    
