
%SLcircBin.m
% 
%  author: Steeve Laquitaine
%    date: 140807
% purpose: bin circular data
%   usage: 
%
%       x = 135:1:225;
%       binStart = 360;   
%       binSize =  3;
%       binEnd = 270;
%       [count,binTag,edgesVal] = SLcircBin(x,binStart,binSize,binEnd)
% 

function [count,binTag,edgesVal] = SLcircBin(x,binStart,binSize,binEnd)

%check
if SLisinteger(binSize)==0
    fprintf('%s /n','binSize must be an integer')
    return
end
if SLisinteger(binStart)==0
    fprintf('%s /n','binStart must be an integer')
    return
end

%linearize egdes space ]-inf + inf[ 
edges = binStart: binSize: binStart + 359;

%convert to circular space in deg
edgesCirc = SLdegLin2Circ(edges);

edgesCirc(edgesCirc==360) = 0;

%stay in range
ccwRange = SLcircCcwVect2Angle(binStart,binEnd);
edgesCircccwDist = SLcircCcwVect2Angle(binStart,edgesCirc);
edgesCircRange = edgesCirc(edgesCircccwDist <= ccwRange);

%sort
edgesCircSort = sort(edgesCircRange);

%1 degree unit
edgesCircSort = round(edgesCircSort);
%rad
%edgesRad = de2r(edges,0);
%xRad = de2r(x,0) - de2r(binStart,1);
[count,binTag] = histc(x,edgesCircSort);
count = count';
binTag = binTag';
edgesVal = edgesCircSort';






