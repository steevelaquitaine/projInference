

%SLtextIsContained.m
%
%  author: steeve laquitaine
%    date: 2015/04/30
% purpose: scan word and find a piece of successive letters in the word
%
%   usage: 
%
%           [match,matchPos] = SLtextIsContained('hello','hel')

function [match,matchPos] = SLtextIsContained(word,piece)

%scan word with piece
for i = 1 : length(word)-length(piece)+1 
    
    %scan window
    scanWnd = i:i+length(piece)-1;
    
    %match
    matchPos(i) = strcmp(word(scanWnd),piece);    
end

%say if match
if sum(matchPos)>0
    match = 1;
else
    match = 0;
end

