%SLdrawModelsPredictionHistwithDataCentered0.m
%
% author: steeve laquitaine
%   date: 141122 updated 150724
%purpose: Draw model predicted estimate distributions
%
%  usage:
%
%       SLdrawModelsPredictionHistwithDataCentered0(Pdata,Bins,Ppred,d,coh,pstd,priorModes,...
%           cond,priorShape,varargin)
%
%Arguments:
%
%       Bins : vector of 1 by Nbins (10:10:360)
%      Pdata : matrix of Nbins (e.g., 36) by Ncond (e.g., 202)
%      Ppred : matrix of Nbins (e.g., 36) by Ncond (e.g., 202)
%          d : series of factor 1 (e.g., motion direction)
%        coh : series of factor 2 (e.g., coherence)
%       pstd : series of factor 3 (e.g., prior std)
% priorModes : series of factor 4 (e.g., prior modes)
%       cond : matrix of Ncond by 3 factors (1 to 3)
% priorShape : prior shape :'vonMisesPrior' or 'bimodalPrior'
%   varargin : {'experiment','vonMisesPrior'} or {'experiment','bimodalPrior'}
%
%
% see SLfitBayesianModel.m


%model distribution (centered at prior mean)
function SLdrawModelsPredictionHistwithDataCentered0(Pdata,Bins,Ppred,d,coh,pstd,priorModes,...
    cond,priorShape,varargin)

%args
varargs = varargin{1};

%case meanData or stdData is not input
if isempty(Pdata)==1
    Pdata = nan(size(Ppred));
end
if isempty(Ppred)==1
    Ppred = nan(size(Pdata));
end

%warning
if isempty(Ppred) && isempty(Pdata)==1
    fprintf('%s /n','(SLdrawModelsPredictionHistwithDataCentered) Both input arguments "Pdata"',...
        ' and "Ppred" are empty. At least one of them should contain data /n')
end

%we plot the bin value in the middle e.g.., bin 0 to 10 is plotted at
%position 5.
middleBinsForPlot = (Bins(2) - Bins(1))/2;
bins = Bins - middleBinsForPlot;

%graphics
colors = [1 0 0];

%priors, coherences and motion directions as signed linear distance to
%prior mean
%case von Mises Prior
%--------------------
if strcmp(priorShape,'vonMisesPrior')
    Thepriors = unique(pstd);
    
    %case bimodal Prior
    %------------------
elseif strcmp(priorShape,'bimodalPrior')
    
    %prior conditions
    priorCond = priorModes(:,2) - priorModes(:,1);
    
    %check that all priors are there
    numPriorCond = size(SLuniqpair(priorModes),1);
    if numel(unique(priorCond)) == numPriorCond
        Thepriors = unique(priorCond);
    end
else
    fprintf('%s \n','(SLdrawModelsPredictionHistwithDataCentered) you need to input prior type /n')
end
ThePriorModes = unique(priorModes);
Thecohs = unique(coh);

%motion directions are centered to the mean of the experimental priors for
%both von Mises and bimodal priors
dlinDisttoPrior = d - nanmean(priorModes,2);

%number of conditions
numPriors = numel(Thepriors);
numCoh = numel(Thecohs);

%set of motion directions as signed linear distance to prior mean
%and associated motion direction in deg
[ThedirlinDisttoPrior,pos] = unique(dlinDisttoPrior);
TheDir = d(pos);
motDir = unique(TheDir);
numDir = numel(motDir);

%case we want to superimpose priors on same axis
%-----------------------------------------------
if sum(strcmp(varargs,'SuperImposedPrior'))==1
    
    %look at how data distribution change as prior strength increases
    %get number of axes and calculate position of plot for each motion
    %direction
    axesPos = 1:1:numel(ThedirlinDisttoPrior)*numPriors;
    dirpos = repmat(1:1:numDir,numel(Thepriors),1);
    dirpos = dirpos';
    dirpos = dirpos(:);
    
    %initialize graphics
    scrsz = get(0,'ScreenSize');
    width = 7;
    F.f2.color={[0.5 0 0],...
        [1 0.2 0],...
        [1 0.6 0],...
        [0.75 0.75 0]};
    
    F.f2.colorPre={[0 0 0],...
        [0.7 0.5 0],...
        [0.8 0.4 0],...
        [0.3 0.3 0]};
    
    %initialize the conditions to plot in each figure (coh)
    condThisCoh = nan(numDir*numPriors,size(cond,2));
    condThisCoh(:,1) = SLreplicateRows(Thepriors,numDir);
    
    %center row axes at prior mean
    shift = ceil(numDir/2 - find(motDir==225));
    DirForShiftedRow = circshift(motDir,[shift 0]);
    condThisCoh(:,3) = repmat(DirForShiftedRow,4,1);
    numAllCond = size(condThisCoh,1);
    
    %put strong prior on top
    [condThisCoh(:,1),IA] = sort(condThisCoh(:,1),'descend');
    condThisCoh(:,3) = condThisCoh(IA,3);
    
    %this coherence
    for j = 1 : numCoh
        
        %figure
        figure('color','w',...
            'Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
        
        %this coherence
        thiscoh = Thecohs(j);
        
        %sync subplots and cond
        condThisCoh(:,2) = repmat(Thecohs(j),numDir*numPriors,1);
        for i = 1 : numAllCond
            
            %time
            %t1 = tic;
            
            %axis
            hs(i) = subplot(numDir,1,dirpos(i));
            
            %condition
            thisCon = SLfindRow(condThisCoh(i,:),cond);
            thisCon = find(thisCon==1);
            
            %case data
            %draw data and prediction distribution
            if thisCon ~= 0
                hold all
                
                %recentering to prior
                %--------------------
                %recenter all data relative to prior
                %calculate distance to prior mean 225 deg
                %sort distance
                bins = SLmakeColumn(bins);
                xCentered = SLvectors2signedAngle(bins,225,'polar');
                [xCenteredSorted,IA] = sort(xCentered,'ascend');
                
                %data
                yDataCentered = Pdata(:,thisCon);
                yDataCenteredSorted = yDataCentered(IA);
                
                %predictions
                yPredCentered = Ppred(:,thisCon);
                yPredCenteredSorted = yPredCentered(IA);
                
                %sort positions and labels for plot
                xTickCentered = SLvectors2signedAngle(bins,225,'polar');
                [xTickCenteredSorted,I] = sort(xTickCentered,'ascend');
                xtickLabel = bins(I);
                
                %plot
                %data
                priorPos = find(cond(thisCon,1)==flipud(Thepriors));
                area(xCenteredSorted,yDataCenteredSorted,'facecolor',...
                    F.f2.color{priorPos},...
                    'edgecolor','none')
                
                %predictions
                %plot(xCenteredSorted,yPredCenteredSorted,'color',...
                %    F.f2.colorPre{priorPos},...
                %    'linesmoothing','on','linewidth',1);
                area(xCenteredSorted,yPredCenteredSorted,'facecolor',...
                    F.f2.color{priorPos},...
                    'edgecolor','none')
                
                %graphics
                %dists=[Pdata(:,thisCon); Ppred(:,thisCon)];
                %if sum(isnan(dists))==numel(dists)
                %end
                %ylim is max of all plotted data for ecah axis for clarity.
                %When we use the max of all axis, it's hard to see.
                PosDataThisAxis = find(SLfindRow(condThisCoh(i,[2 3]),...
                    cond(:,[2 3]))==1);
                maxPlot(i) = max(max([Pdata(:,PosDataThisAxis) ...
                    Ppred(:,PosDataThisAxis)]));
                
                xlim([-180 +180])
                %ylim([0 maxPlot(i)])
                set(gca,'ytick',[0 fix(maxPlot(i)*10000)/10000],...
                    'ytickLabel',[0 fix(maxPlot(i)*10000)/10000],...
                    'FontName','Helvetica',...
                    'FontWeight','light',...
                    'FontAngle','oblique')
                
                %motion direction
                motDirCentered = SLvectors2signedAngle(cond(thisCon,3),225,'polar');
                plot([motDirCentered motDirCentered],[0 0.5*maxPlot(i)],...
                    '-','color',[0.5 0.5 0.5],'linewidth',2)
                
                %prior mode
                priorModeCentered = SLvectors2signedAngle(ThePriorModes,225,'polar');
                plot([priorModeCentered priorModeCentered],[0 0.5*maxPlot(i)],...
                    '-','color',[0.3 0.6 0.9],'linewidth',2);
                
                %remove axis
                axis off
                set(gca,'ytick',0,'yticklabel',0)
            else
                axis off
            end
            box off
            drawnow
            %toc(t1)
            %thetime(i) = toc(t1);
            
            %more graphics
            %-------------
            %title top axes
            if dirpos(i)==1
                myprior = Thepriors(dirpos(i));
                
                title(['(',num2str(thiscoh*100),'% coh -',...
                    num2str(myprior),' deg std prior)'],...
                    'fontsize',11)
            end
            
            %no x-tick
            set(gca,'xtick',[],'xtickLabel',[],'fontsize',8,...
                'FontName','Helvetica',...
                'FontWeight','light',...
                'FontAngle','oblique')
            
            %ylabel left axes
            if sum(dirpos(i)==dirpos(1:numDir))==1
                ylabel('Probability','fontsize',11)
            end
            
            %xlabel bottom axes
            if sum(dirpos(i)==numDir)==1
                axis on
                xlabel('Estimated directions (deg)','fontsize',11)
                
                %xmin and xmax
                xlim([-180 +180])
                
                %labels
                set(gca,'xtick',xTickCenteredSorted(3:8:end),...
                    'xticklabel',xtickLabel(3:8:end),...
                    'FontName','Helvetica',...
                    'FontWeight','light',...
                    'FontAngle','oblique')
            end
        end
        set(hs,'ylim',[0 max(maxPlot)])
        SLremoveDeadSpace
    end
    
else
    
    %case prior conditions are sorted in column axes (default)
    %------------------------------------------------------------
    %look at how data distribution change as prior strength increases
    %get number of axes and calculate position of plot for each motion
    %direction
    axesPos = 1:1:numel(ThedirlinDisttoPrior)*numPriors;
    dirpos = [];
    for i = 1 : numel(Thepriors)
        dirpos = [dirpos; i:numel(Thepriors):numel(axesPos)];
    end
    dirpos = dirpos';
    dirpos = dirpos(:);
    
    %initialize graphics
    scrsz = get(0,'ScreenSize');
    width = 7;
    
    F.f2.color={[0.5 0 0],...
        [1 0.2 0],...
        [1 0.6 0],...
        [0.75 0.75 0]};
    
    F.f2.colorPre={[0 0 0],...
        [0.6 0.1 0],...
        [0.8 0.4 0],...
        [0.3 0.3 0]};
    
    %initialize the conditions to plot in each figure (coh)
    condThisCoh = nan(numDir*numPriors,size(cond,2));
    condThisCoh(:,1) = SLreplicateRows(Thepriors,numDir);
    
    %center row axes at prior mean
    shift = ceil(numDir/2 - find(motDir==225));
    DirForShiftedRow = circshift(motDir,[shift 0]);
    condThisCoh(:,3) = repmat(DirForShiftedRow,numPriors,1);
    numAllCond = size(condThisCoh,1);
    %condThisCoh(:,3) = repmat(motDir,4,1);
    %numAllCond = size(condThisCoh,1);
    
    %histogram positions
    figPos = [402 60 434 768; 688 61 206 768; 894 61 206 768];
    
    %this coherence
    for j = 1 : numCoh
        
        %figure
        figure('color','w',...
            'Position',figPos(j,:))
        
        %this coherence
        thiscoh = Thecohs(j);
        
        %sync subplots and cond
        condThisCoh(:,2) = repmat(Thecohs(j),numDir*numPriors,1);
        for i = 1 : numAllCond
            
            %time
            t1 = tic;
            
            %axis
            hs(i) = subplot(numDir,numPriors,dirpos(i));
            
            %condition
            thisCon = SLfindRow(condThisCoh(i,:),cond);
            thisCon = find(thisCon==1);
            
            %case data
            %draw data and prediction distribution
            if thisCon ~= 0
                hold all
                
                %recentering to prior
                %--------------------
                %recenter all data relative to prior
                %calculate distance to prior mean 225 deg
                %sort distance
                xCentered = SLvectors2signedAngle(bins,225,'polar');
                [xCenteredSorted,IA] = sort(xCentered,'ascend');
                
                %data
                yDataCentered = Pdata(:,thisCon);
                yDataCenteredSorted = yDataCentered(IA);
                
                %predictions
                yPredCentered = Ppred(:,thisCon);
                yPredCenteredSorted = yPredCentered(IA);
                
                %sort positions and labels for plot
                xTickCentered = SLvectors2signedAngle(bins,225,'polar');
                [xTickCenteredSorted,I] = sort(xTickCentered,'ascend');
                xtickLabel = bins(I);
                
                %plot
                %data
                priorPos = find(cond(thisCon,1)==flipud(Thepriors));
                area(xCenteredSorted,yDataCenteredSorted,'facecolor',...
                    F.f2.color{priorPos},...
                    'edgecolor','none')
                
                %predictions
                plot(xCenteredSorted,yPredCenteredSorted,'color',...
                    F.f2.colorPre{priorPos},...
                    'linesmoothing','on','linewidth',1);
                
                %graphics
                %dists=[Pdata(:,thisCon); Ppred(:,thisCon)];
                %if sum(isnan(dists))==numel(dists)
                %end
                %ylim is max of all plotted data for ecah axis for clarity.
                %When we use the max of all axis, it's hard to see.
                maxPlot(i) = max([yDataCenteredSorted; yPredCenteredSorted]);
                xlim([-180 +180])
                ylim([0 maxPlot(i)])
                set(gca,'ytick',[0 fix(maxPlot(i)*1e6)/1e6],...
                    'ytickLabel',[0 fix(maxPlot(i)*1e6)/1e6],...
                    'FontName','Helvetica',...
                    'FontWeight','light',...
                    'FontAngle','oblique',...
                    'linewidth',0.5)
                
                %motion direction marker
                motDirCentered = SLvectors2signedAngle(cond(thisCon,3),225,'polar');
                plot([motDirCentered motDirCentered],[0 0.25*maxPlot(i)],...
                    '-','color',[0.8 0.8 0.8],'linewidth',2)
                
                %prior modes marker
                %prior modes in this condition
                if strcmp(priorShape,'bimodalPrior')
                    ThePriorModes = SLuniqpair(priorModes(priorCond==cond(thisCon,1),:));
                elseif strcmp(priorShape,'vonMisesPrior')
                end
                priorModeCentered = SLvectors2signedAngle(ThePriorModes,225,'polar');
                clear ijk
                for ijk = 1 : numel(priorModeCentered)
                    plot([priorModeCentered(ijk) priorModeCentered(ijk)],[0 0.25*maxPlot(i)],...
                        '-','color',[0.3 0.6 0.9],'linewidth',2);
                end
                %remove axis
                axis off
                set(gca,'ytick',0,'yticklabel',0,...
                    'FontName','Helvetica',...
                    'FontWeight','light',...
                    'FontAngle','oblique')
            else
                axis off
            end
            box off
            drawnow
            %toc(t1)
            %thetime(i) = toc(t1);
            
            %more graphics
            %-------------
            %title top axes
            if sum(dirpos(i)==1:numPriors)==1
                myprior = Thepriors(dirpos(i));
                
                title(['(',num2str(thiscoh*100),'% coh -',...
                    num2str(myprior),' deg std prior)'],...
                    'fontsize',11)
            end
            
            %no x-tick
            set(gca,'xtick',[],'xtickLabel',[],'fontsize',10)
            
            %ylabel left axes
            if sum(dirpos(i)==dirpos(1:numDir))==1
                ylabel('Probability','fontsize',11)
            end
            
            %xlabel bottom axes
            bottomAxes = numAllCond-numPriors+1:numAllCond;
            if sum(dirpos(i)==bottomAxes)==1
                axis on
                xlabel('Estimated directions (deg)','fontsize',11)
                
                %xmin and xmax
                xlim([-180 +180])
                
                %labels
                set(gca,'xtick',25:40:355,'xtickLabel',25:40:355,...
                    'FontName','Helvetica',...
                    'FontWeight','light',...
                    'FontAngle','oblique')
                
                set(gca,'xtick',xTickCenteredSorted(3:8:end),...
                    'xticklabel',xtickLabel(3:8:end))
                
                %                 %case von Mises prior
                %                 %--------------------
                %                 if strcmp(priorShape,'vonMisesPrior')
                %                     set(gca,'xtick',xTickCenteredSorted(3:8:end),...
                %                         'xticklabel',xtickLabel(3:8:end))
                %
                %                     %case bimodal prior
                %                     %------------------
                %                 elseif strcmp(priorShape,'bimodalPrior')
                %
                %                     %display prior modes
                %                     [~,pos] = intersect(xtickLabel,unique(priorModes));
                %                     set(gca,'xtick',xTickCenteredSorted(pos),...
                %                         'xticklabel',unique(priorModes))
                %
                %                     %Warning
                %                 else
                %                     fprintf('%s /n','(SLdrawModelsPredictionCentered) you need to input prior type')
                %                 end
            end
        end
        %set(hs,'ylim',[0 max(maxPlot)])
        %clear up
        SLremoveDeadSpace
    end
end