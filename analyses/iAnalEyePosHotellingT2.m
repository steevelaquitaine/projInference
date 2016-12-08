
%iAnalEyePosHotellingT2.m
%
%
%author : steeve laquitaine
%purpose: run hotelling T2 test to test whether a subject's eye positions are
%         significantly different (p<0.05) between prior conditions
%
%usage :
% 
%     irootpath = '~/Dropbox/myDropBox/Codes/';
%     subject = 'sub01';
%     iAnalEyePosHotellingT2(subject,irootpath)


function iAnalEyePosHotellingT2(subject,irootpath)

%project paths
mungedatapath = [irootpath 'projInference/data/eye/'];
subjectpath = [mungedatapath subject];
cd(subjectpath)
load eyexy.mat
load priorstd.mat

%remove missing data
[missing,col] = find(isnan(eyexy(:,1)));
eyexy(missing,:) = [];
priorstd(missing,:) = [];

%sort eye positions by prior and do pairwise hotelling T^22 testhote
%convert x,y to angles in radians and distance to fixation (length)
pu = unique(priorstd);
for i = 1 : length(pu)
   eyeByPrior{i} = eyexy(priorstd==pu(i),:);    
   [eyePosRadP{i},~,eyeDistToFixP{i}] = SLcart2polar(eyeByPrior{i});
end

%hotelling t2 test for each pair of prior
close all
pairsOfPrior = nchoosek(pu,2); 
pairsOfPriorIx = nchoosek(1:length(pu),2);
for i = 1 : length(pairsOfPrior)    
    eyePosRadP1 = eyePosRadP{pairsOfPriorIx(i,1)}; 
    eyePosRadP2 = eyePosRadP{pairsOfPriorIx(i,2)};
    eyeDistToFixP1 = eyeDistToFixP{pairsOfPriorIx(i,1)}; 
    eyeDistToFixP2 = eyeDistToFixP{pairsOfPriorIx(i,2)};   
    %hotelling t2 for pairs of priors         
    p(i) = hotelling(eyePosRadP1,eyePosRadP2,eyeDistToFixP1,eyeDistToFixP2,1);
    title({sprintf('%s %i %s %i','Priors',pairsOfPrior(i,1),'vs.',pairsOfPrior(i,2)),...
        sprintf('%s %i','Hotelling t2 p-value = ',p(i))})        
end

%get one clear figure
figs = sort(SLgetFigures);
figure(max(figs+1));
set(gcf,'color','w')
for i = 1 : max(figs)
    figure(i);
    oldAxis = gca;
    figure(max(figs+1));       
    newAxis = subplot(1,max(figs),i);
    slgraphCopyFigAxtoNewfigAx(oldAxis,figure(max(figs+1)),newAxis)
    if i == 1
        hold all; fix=plot(0,0,'.k');
        xlabel('Horizontal position (deg)')
        ylabel('Vertical position (deg)')        
        legend(fix,'Fixation point at (0,0)')
    end    
end

SLConventionUp
close(figs)
