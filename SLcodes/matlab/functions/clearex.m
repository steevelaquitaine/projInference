
 %usage: clearex('fileName1','fileName2')
%Author: Arnaud Laurent
  %Date: Feb 15th 2008

%Note: It is possible to use wildcard (e.g. 'a*')


function clearex(varargin)
%% This function clear all workspace
% except for one or more selected variable

a = evalin('base','who');

var = cell(size(varargin));

for i=1:nargin
    var{i}=varargin{i};
end

assignin('base','ClEaRGsJioU',var);
var = evalin('base','who(ClEaRGsJioU{:})');

clearvar = a(~ismember(a,var));
assignin('base','ClEaRGsJioU',clearvar);

evalin('base','clear(ClEaRGsJioU{:},''ClEaRGsJioU'')')

