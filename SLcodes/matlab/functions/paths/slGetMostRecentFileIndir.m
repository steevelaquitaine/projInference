
%slGetMostRecentFileIndir.m
%
%
% author: steeve laquitaine
%   date: 160121
%purpose: get most recent .mat file in directory
%
%  usage:
%
%       lastfile = slGetMostRecentFileIndir
%
%    ref: https://www.mathworks.com/matlabcentral/newsreader/view_thread/270611


function lastfile = slGetMostRecentFileIndir

d = dir('*.mat');
[dx,dx] = sort([d.datenum],'descend');
lastfile = d(dx==1).name;