%SLuniqpair.m
%
%       $Id: SLuniqpair.m 750 2014-05-15 12:23:46Z steeve $
%     usage: analyses
%        by: steeve laquitaine
%      date: 140513
% copyright: (c) 2014 steeve laquitaine
%   purpose: find unique rows combinations in a matrix
%             
%     usage: 
%        pairs=[1 1 1 1 1 1 1 1; 1 1 1 1 2 2 2 2; 1 2 3 3 1 2 3 3]';
%        uniquePair=SLuniqpair(pairs);
%     
%Description:
%   uniquePair are the unique combination
%   b are the rows indices for each combinations in pairs
%   c are the rows indices for each combinations in uniquePair

function [uniquePair,b,c] = SLuniqpair(pairs)

%[uniquePair,b,c]=unique(sort(pairs,2), 'rows');
[uniquePair,b,c] = unique(pairs,'rows');
