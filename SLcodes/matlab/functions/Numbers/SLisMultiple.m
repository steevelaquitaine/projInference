

% author: steeve laquitaine
%   date: 150415
%purpose: find if a number is a multiple of another number 

%  usage:
%
%           yes = SLisMultiple(x2,x1)

function yes = SLisMultiple(x2,x1)

yes = (mod(x2,x1)==0);
