
%memo05.m
%SID: 322
%
%
%naive in psychophysics
%first mglSetSID('SID') then run line by line
%
%
%--------------------------- notes -------------------------------------
%First tried with con[.1 .55 1]
%--> 100% contrasts shows sub understands task
%--> 55 and 10% contrasts shows sub bias towards prior mean and variab
%    increases
%Search intermediate contrast to get intermediate bias and variability
%during session 1
%
%Adjusted intermediate contrast to get intermediate bias and variability
%at con[0.1 - 0.156 -1] starting from session 2
%-----------------------------------------------------------------------
%
%searching good contrast
%-----------------------
%1-151208
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.55 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.3 -1]
% 
%
%final experiment
%----------------
%1-151208
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%[0.1 - 0.156 -1]
%
%2-151210
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%[0.1 - 0.156 -1]
%
%3-151211
%sltaskLoc('displayName=VPixx','computeStim','p=40');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%
%4-160108
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%
%5-160129
%sltaskLoc('displayName=VPixx','computeStim','p=40');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%
%6-160129
%sltaskLoc('displayName=VPixx','computeStim','p=10');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%[0.1 - 0.156 -1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%[0.1 - 0.156 -1]
