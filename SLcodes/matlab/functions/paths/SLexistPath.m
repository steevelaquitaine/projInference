

% author: steeve laquitaine
%   date: 15/04/01
%purpose: check if path exists
%
%usage:
%
%       SLexistPath('/Users/steeve_laquitaine/Dropbox/')


function output = SLexistPath(path)

output = exist(path)>0;

sprintf(['I found your path: ',path])