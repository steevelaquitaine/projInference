

%SLLoadTextData.m
%
% author: Steeve Laquitaine
%   date: 141229
%purpose: Load text data
% 
%  usage:  SLLoadTextData('Tutorial')
%
%options (varargin)


function SLLoadTextData


%case tutorial
if nargin==1 && strcmp(x,'Tutorial')
    
    %open tutorial data file
    o.fid = fopen('irisdata','r');
    
    %read data
    o.datatxt = textscan(o.fid,'%s'); 
    o.datatxt = o.datatxt{:};
    
    %nb rows
    o.nlines = size(o.datatxt,1);
    
    %columns delimits
    for i = 1 : o.nlines
        cdelim = find(o.datatxt{i}==',');
        cstart{i} = [1 cdelim + 1];
        cends{i} =  [find(o.datatxt{i}==',')-1 size(o.datatxt{i},2);];
    end    
    
    o.ncol = numel(cstart{1});
    
    %isolate all values in a cell
    for i = 1 : o.nlines
        for j = 1 : o.ncol
            o.data{i,j} = o.datatxt{i}(cstart{i}(j) : 1: cends{i}(j));
        end
    end
end