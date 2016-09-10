
%SLisXLabel.m
%
% author: Steeve laquitaine
%   date: 140712
%purpose: say if an axis has an Xlabel (1) or not (0)
%  usage: SLisXLabel(axhdle1, axhdle2)


function output = SLisXLabel(axhdle)

%call for help
if ieNotDefined('axhdle')
    help SLisXLabel
  return    
end

%There's an Xlabel if the height is non zero
if  iscell(get(axhdle,'xlabel'))
    position = cell2mat(get(cell2mat(get(axhdle,'xlabel')),'extent'));
    output = position(:,4)>0;
else
    position = get(get(axhdle,'xlabel'),'extent');
    output = position(:,4)>0;
end