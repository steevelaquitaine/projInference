


%slfMRIplotVoxTuningVsChance.m
%
% author: steeve laquitaine
%purpose: plot permuted null vectors calculated for a tuning function with 
%         Georgopoulos popultion vector analysis (by shuffling responses 
%         and associated variable) and observed vector
%
%  usage: 
%
%       slfMRIplotVoxTuningVsChance(nullVecsDir,obsVeclen,obsVec,meanVlenBydir,VlenCIbyDir,vecdir)       
%
% 


function slfMRIplotVoxTuningVsChance(nullVecsDir,obsVeclen,obsVec,meanVlenBydir,VlenCIbyDir,vecdir)

%null (grey) and observed (red) vectors
subplot(1,2,1)
maxplot = max([obsVeclen; meanVlenBydir]);
sldrawCircErrorBar(vecdir,meanVlenBydir,VlenCIbyDir,maxplot);
hold on; h = polar(round(obsVec.deg.mean),obsVeclen,'ro');
set(h,'markersize',8)

%circular pdf
subplot(1,2,2)
vecdirPdf = hist(nullVecsDir,1:1:360);
sldrawVecPdf(1:1:360,vecdirPdf)