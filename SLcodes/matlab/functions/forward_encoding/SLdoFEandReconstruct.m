%SLdoFEandReconstruct
% 
% Author: Steeve Laquitaine
%   date: 140411
%purpose:
%         run forward encoding (FE)
%         reconstruct parameter from trained weights-->test channels
%         outputs and known channels' output.
%
%  usage: 
%
%       [reconstParam] = SLdoFEandReconstruct(BOLDmnTrain,BOLDmnTest,...
%       paramValsTrain,paramValsTest,paramSelectivity,OptionDisp)

         
function [reconstParam,Condition,R,P] = SLdoFEandReconstruct(BOLDmnTrain,BOLDmnTest,paramValsTrain,paramValsTest,paramSelectivity,OptionDisp)

%do forward encoding
[Ctrkn,~,~,~,Wmk]=SLdoForwardEncoding(BOLDmnTrain,...,
    6,...
    paramValsTrain,...
    OptionDisp);

%display option
if strcmp(OptionDisp,'display=on')==1
    h=get(gcf,'children');
    set(h(3),'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
end

%reconstruction of test parameters (e.g., motion directions). We use the
%trained weights to generate the channels output in test data for each
%trial. Then the reconstructed parameter is the parameter that have knwonfor each trial we find the known channel output that 
%correlates best with the test trial channel output.
[reconstParam,~,~,Condition]=reconstructDir(BOLDmnTest,Wmk,Ctrkn,paramValsTrain);

%display option
if strcmp(OptionDisp,'display=on')==1
    figure; drawCircStat(reconstParam',paramValsTest'); close;
    plot([0 360],[225 225],'--b')
    title(['Performance (R=',num2str(R(2,1)),') (p=',...
    num2str(P(2,1)),') (Pearson corr)'],'fontsize',14)
end

%get decoding/reconstruction performances
[R,P]=corrcoef(reconstParam,paramValsTest);

%Warnings
fprintf('%12s \n',['IMPORTANT WARNING: be careful, weight matrix may be',...
    'ill-conditioned (Condition=',num2str(Condition),')'])
fprintf('%12s \n','Decoding of parameters...done')
fprintf('%12s \n',['Performance (R=',num2str(R(2,1)),') (p=',...
    num2str(P(2,1)),') (Pearson corr)'])

%Visualize the data that we model with forward encoding. What we want to 
%see is the average BOLD response of each voxel (m row) for each motion 
%direction displayed. We need to average voxel response for repetition of 
%the same motion direction.
%training data
paramsValuniq=unique(paramValsTrain);
numparamsValuniq=numel(paramsValuniq);
mall=nan(size(BOLDmnTrain,1),numparamsValuniq);
sall=nan(size(BOLDmnTrain,1),numparamsValuniq);
for j=1:size(BOLDmnTrain,1)
    [m,s]=makeStat(BOLDmnTrain(j,:)',paramValsTrain);
    mall(j,:)=m;
    sall(j,:)=s;
end

if strcmp(OptionDisp,'display=on')==1
    figure('color','w');
    subplot(121)
    imagesc(mall)
    colorbar
    set(gca,'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
    set(gca,'xtick',1:numparamsValuniq,'xticklabel',paramsValuniq)
    ylabel({'Voxels sorted by','direction selectivity (degrees)'},'fontsize',14)
    title('Raw BOLD (Training)','fontsize',14)
    colormap('hot')
    axis square
    xlabel('Motion direction (degrees)','fontsize',14)
end

%test data
paramsValuniq=unique(paramValsTest);
numparamsValuniq=numel(paramsValuniq);
mall=nan(size(BOLDmnTest,1),numparamsValuniq);
sall=nan(size(BOLDmnTest,1),numparamsValuniq);
for j=1:size(BOLDmnTest,1)
    [m,s]=makeStat(BOLDmnTest(j,:)',paramValsTest);
    mall(j,:)=m;
    sall(j,:)=s;
end

if strcmp(OptionDisp,'display=on')==1
    subplot(122)
    imagesc(mall)
    colorbar
    set(gca,'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
    set(gca,'xtick',1:numparamsValuniq,'xticklabel',paramsValuniq)
    ylabel({'Voxels sorted by','direction selectivity (degrees)'},'fontsize',14)
    title('raw BOLD (Test)','fontsize',14)
    colormap('hot')
    axis square
    xlabel('Motion direction (degrees)','fontsize',14)
end

%We can also look at raw and filtered data.
%raw
if strcmp(OptionDisp,'display=on')==1
    figure('color','w');
    subplot(121)
    imagesc(BOLDmnTrain)
    colorbar
    set(gca,'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
    set(gca,'xtick',1:30:numel(paramValsTrain),'xticklabel',paramValsTrain(1:30:end))
    ylabel({'Voxels sorted by','direction selectivity (degrees)'},'fontsize',14)
    title('raw BOLD (Training)','fontsize',14)
    colormap('hot')
    axis square
    
    subplot(122)
    imagesc(BOLDmnTest);
    colorbar
    SLpositionFigure(gcf,1,1)
    set(gca,'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
    set(gca,'xtick',1:30:numel(paramValsTest),'xticklabel',paramValsTest(1:30:end))
    title('raw BOLD (Test)','fontsize',14)
    colormap('hot')
    axis square
end
% %smoothed with a 2-d filter. We only smooth across columns (direction) not
% %across rows (voxels). We need to denoise voxel activity across directions
% subplot(223)
% %myfilter=fspecial('gaussian',[1 10], 0.5);
% myfilter=fspecial('average',[1 7]);
% BOLDmnTrainfilt=imfilter(BOLDmnTrain,myfilter,'replicate');
% imagesc(BOLDmnTrainfilt)
% set(gca,'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
% set(gca,'xtick',1:30:numel(paramValsTrain),'xticklabel',paramValsTrain(1:30:end))
% colorbar
% title('Filtered BOLD','fontsize',14)
% ylabel({'Voxels sorted by','direction selectivity (degrees)'},'fontsize',14)
% xlabel('Motion direction (degrees)','fontsize',14)
% colormap('hot')
% axis square
% 
% subplot(224)
% BOLDmnTestfilt=imfilter(BOLDmnTest,myfilter,'replicate');
% imagesc(BOLDmnTestfilt)
% colorbar
% set(gca,'ytick',1:10:numel(paramSelectivity),'yticklabel',paramSelectivity(1:10:end))
% set(gca,'xtick',1:30:numel(paramValsTest),'xticklabel',paramValsTest(1:30:end))
% title('Filtered BOLD (Test)')
% xlabel('Motion direction (degrees)','fontsize',14)
% colormap('hot')
% axis square









