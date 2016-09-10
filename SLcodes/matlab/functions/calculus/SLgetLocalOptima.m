
%SLgetLocalOptima.m
%
% author: steeve laquitaine
%   date: 150631
%purpose: get a function local optima (where slopes changes)
%  usage:
%
%       pos = SLgetLocalOptima(x)

function pos = SLgetLocalOptima(x,type)

signx = sign(x);
pos = find(signx(2:end) + signx(1:end-1) == 0); %rows where slope changes

if strcmp(type,'minima');
    
    minima = find(x(pos) - x(pos+1) < 0);
    pos = pos(minima);
    
elseif strcmp(type, 'maxima');
    
  maxima = find(x(pos) - x(pos+1) > 0);
    pos = pos(maxima);
end
