%slCalcMeanByRowCondEachCol.m
%
%
% author: steeve laquitaine
%purpose: calculate the means over conditions labelled in the rows of an 
%         m by n matrix
%       
%e.g., if 1 (up), 2 (stable) and 3(down) are three conditions we get:
%
%       rowCond        data
%       
%          1        a  i  q  y
%          1        b  j  .
%          2        c  k  .
%          2        d  .
%          3        e  .
%          3        f
%          1        g
%          1        h
%   
%
%usage : 
%
%       data = randi(8,4);
%       rowCond = [1 1 2 2 3 3 1 1];  %contains conditions
%       meanbyCond = slCalcMeanByRowCondEachCol(data,rowCond);
%
%
%    rowCond : vector of labelled conditions
%       data : matrix of rows by columns

function meanbyCond = slCalcMeanByRowCondEachCol(data,rowCond)

nCond = length(unique(rowCond));
data = data';
[row_loc, col_loc] = ndgrid(rowCond,1:nCol);
meanbyCond = accumarray([row_loc(:) col_loc(:)],data(:),[nCond nperm],@mean);
