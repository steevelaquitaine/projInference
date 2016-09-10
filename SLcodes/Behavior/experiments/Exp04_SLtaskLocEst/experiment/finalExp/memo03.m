
%memo03.m
%SID: 320
%
%
%naive in pychophysics
%first mglSetSID('SID') then run line by line
%
%--------------------------- notes -------------------------------------
%First tried with con[.1 .12 1] 
%--> 100% contrasts shows sub understands task
%--> 10 and 12% contrasts shows sub bias towards prior mean and variab
%    increases
%
%Adjusted intermediate contrast to get intermediate bias and variability
%con[.1 .156 1] starting session 2
%-----------------------------------------------------------------------
%
%piloting contrast
%----------------
%1-151207
%sltaskLoc('displayName=VPixx','computeStim','p=80');con[.1 .12 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');con[.1 .12 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');con[.1 .12 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');con[.1 .12 1]
%
%
%final Experiment
%----------------
%1-151208 (very good keep those parameters)
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%
%2-151209
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%
%3-151209
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%
%4-151210
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%
%5-151210
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
  