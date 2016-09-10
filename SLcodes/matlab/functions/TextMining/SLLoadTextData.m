

%SLLoadTextData.m
%
% author: Steeve Laquitaine
%   date: 141229
%purpose: Load text data
%
%  usage:
%
%       o = SLLoadTextData('Tutorial')
%
%       o = SLLoadTextData('Users/steeve/Desktop/myDatafile.txt')
%
%
%Description:
%
%   x must be a text file (.txt) where columns either contain only
%   strings or numbers

function o = SLLoadTextData(x)


%case tutorial
if strcmp(x,'Tutorial')
    
    %open tutorial data file
    o.fid = fopen('irisdata','r');
    
else
    %open input text data file
    o.fid = fopen(x,'r');
end

%read data
o.datatxt = textscan(o.fid,'%s');
o.datatxt = o.datatxt{:};
fclose(o.fid);

%nb rows (instances)
o.nlines = size(o.datatxt,1);

%columns delimits
for i = 1 : o.nlines
    cdelim = find(o.datatxt{i}==',');
    cstart{i} = [1 cdelim + 1];
    cends{i} =  [find(o.datatxt{i}==',')-1 size(o.datatxt{i},2);];
end

%nb of cols (dimensions)
o.ncol = numel(cstart{1});

%isolate all values in a cell
for i = 1 : o.nlines
    for j = 1 : o.ncol
        o.data{i,j} = o.datatxt{i}(cstart{i}(j) : 1: cends{i}(j));
    end
end

%Isolate quantitative and qualitative dimensions of the data set.
for i = 1 : o.nlines
    for j = 1 : o.ncol
        num0{i,j} = str2num(o.data{i,j});
    end
end
for j = 1 : o.ncol
    strPos(j) = isempty(num0{1,j});
    numPos(j) = ~isempty(num0{1,j});
end
o.qualData = o.data(:,strPos);
o.quanData = cell2mat(num0(:,numPos));
