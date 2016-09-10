
%SLdegLin2Circ.m
% 
%  author: Steeve Laquitaine
%    date: 140807
% purpose: convert angle on the linear space ]-inf +inf[ in deg in the
%          circular space [1 360]
%
%   usage: 
%
%       degCirc = SLdegLin2Circ([1 270 370 450]);

function degCirc = SLdegLin2Circ(degLin)

lin2rad = de2r(degLin,0);
degCirc = ra2d(lin2rad);