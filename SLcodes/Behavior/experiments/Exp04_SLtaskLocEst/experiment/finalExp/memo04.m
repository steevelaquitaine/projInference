
%memo04.m
%SID: 321
%
%
%naive in psychophysics
%first mglSetSID('SID') then run line by line
%
%
%--------------------------- notes -------------------------------------
%First tried with con[.1 .12 1] 
%--> 100% contrasts shows sub understands task
%--> 10 and 12% contrasts shows sub bias towards prior mean and variab
%    increases but by same amount
%
%Adjusted intermediate contrast to get intermediate bias and variability
%con[.1 .12 1] during session 2
%
%fixed intermediate contrast --> intermediate bias
%con[.1 .156 1] starting session 3
%-----------------------------------------------------------------------
%
%piloting
%--------
%1-151207
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .12 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .12 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .12 1]
%
%2-151208
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .2 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .15 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .175 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .1375 1]
%
%
%Final experiment
%----------------
%1-151209
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%
%2-151209
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%
%3-151211
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=10');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%
%4-160108
%sltaskLoc('displayName=VPixx','computeStim','p=40');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1]
%sltaskLoc('displayName=VPixx','computeStim','p=20');%con[.1 .156 1] %subject is clearly falling asleep -- discard !!!!
%sltaskLoc('displayName=VPixx','computeStim','p=80');%con[.1 .156 1] %good 
