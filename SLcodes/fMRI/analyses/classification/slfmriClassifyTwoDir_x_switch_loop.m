
%slfmriClassifyTwoDir_x_switch_loop
%
%
% author: steeve laquitaine
%purpose: loop over ROIs and classify two motion directions from voxel 
%         patterns sorted by switching for 6% coherence
%
%usage:
%
%       slfmriInitClassifAnalysisTaskDotDirfMRI05_loop
%       [accuracy,ste,dataClassifsw1,dataClassifsw2,nTrials] = slfmriClassifyTwoDir_x_switch_loop(o);


function [accuracy,ste,dataClassifsw1,dataClassifsw2,nTrials] = slfmriClassifyTwoDir_x_switch_loop(o)

%loop over rois
for i = 1 : length(o.myROInameAll) 
    
    %set this roi
    o.myROIname = {o.myROInameAll{i}};            
    
    %classify
    [dataClassifsw1{i},dataClassifsw2{i}] = slfmriClassifyTwoDir_x_switch(o);
        
    %collect data
    accuracy(i,1) = dataClassifsw1{i}.c.myClasf.raw.fullSg.correct;
    ste(i,1) = dataClassifsw1{i}.c.myClasf.raw.fullSg.correctSTE;
    nTrials.sw1dir1 = length(dataClassifsw1{i}.c.myClasf.raw.fullSg.classifierOut{1});
    nTrials.sw1dir2 = length(dataClassifsw1{i}.c.myClasf.raw.fullSg.classifierOut{2});
    
    accuracy(i,2) = dataClassifsw2{i}.c.myClasf.raw.fullSg.correct;
    ste(i,2) = dataClassifsw2{i}.c.myClasf.raw.fullSg.correctSTE;
    nTrials.sw2dir1 = length(dataClassifsw2{i}.c.myClasf.raw.fullSg.classifierOut{1});
    nTrials.sw2dir2 = length(dataClassifsw2{i}.c.myClasf.raw.fullSg.classifierOut{2});
end

figure('color','w')

subplot(1,2,1)
bar(accuracy(:,1),'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(accuracy(:,1)),accuracy(:,1),'yError',ste(:,1),'Symbol=.','MarkerSize=1')
hline(0.5,'--r')
set(gca,'xtick',1:length(accuracy(:,1)),'xticklabel',o.myROInameAll)
rotateXLabels(gca(),60)
axis tight
ylim([0 1])
xlim([0 length(accuracy)+1])
box off
title('Decode directions 15 and 85 switching to prior (cv ste)')

subplot(1,2,2)
bar(accuracy(:,2),'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(accuracy(:,2)),accuracy(:,2),'yError',ste(:,2),'Symbol=.','MarkerSize=1')
hline(0.5,'--r')
set(gca,'xtick',1:length(accuracy(:,2)),'xticklabel',o.myROInameAll)
rotateXLabels(gca(),60)
axis tight
ylim([0 1])
xlim([0 length(accuracy)+1])
box off
title('Decode directions 15 and 85 switching to motion direction (cv ste)')