
%slStatAncova.m
%
%
% author: steeve laquitaine
%purpose: do analysis of covariance
%
%Input:
%
%       x: displayed dir
%       y: estimated dir
%   group: "strong prior" vs "weak prior"
%
%
%OUTPUT: 
%
%   Stat : interaction p-value in Anova table specifies if slopes differ.
%
%note:
%%Important issue! I can't find how to apply ancova on circular data.
%Solution: Ancova can be performed with the vector averaged 
%estimated dirs and the displayed dirs.

function [Stat] = slStatAncova(x,y,group)

%check enough group
if length(unique(group))<2
    Stat = [];
    sprintf('(statAncova) Not enough groups for ANCOVA')
    return
else
    [h,atab,ctab,stats] = aoctool(x,y,group,0.05,'','','','off');
    Stat = {'h','anovatab','coeftab','stats'};
    Stat(2,:) = {h,atab,ctab,stats};
end
