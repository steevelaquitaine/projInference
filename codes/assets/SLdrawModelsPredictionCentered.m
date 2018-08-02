
%SLdrawModelsPredictionCentered.m
%
% author: steeve laquitaine
%   date: 141122 updated 150710
%purpose: Draw models predictions mean and std centered with data
%
%usage:
%
%       meanData : N by 1 vector of data
%       stdData  : N by 1 vector of data std
%       meanPred : N by 1 vector of pred
%       stdPred  : N by 1 vector of pred std
%       dataCond : N by 3 matrix of conditions for each value of
%       priorShape : 'vonMisesPrior' or 'bimodalPrior'
%       priorModes : if priorShape is 'vonMisesPrior' --> one scalar 
%                  : vector of two values 
%
%varagin
%
%       'yCentered' :  Center y axis
%                 
%       SLdrawModelsPredictionCentered(meanData,stdData,meanPred,errorPred,dataCond,...
%           priorModes,priorShape,'yCentered')
%
%Description:
%
%       see SLfitBayesianModel

function SLdrawModelsPredictionCentered(meanData,errorOfmeanData,stdData,...
    errorOfstdData,meanPred,errorOfmeanPred,stdPred,errorOfstdPred,dataCond,...
    priorModes,priorShape,varargin)

%-----------------------------------------
%checks
%-----------------------------------------
%check varargin
if isempty(varargin)
    varargin = {'empty'};
end

%make sure column vector
meanPred = SLmakeColumn(meanPred);
meanData = SLmakeColumn(meanData);
stdData = SLmakeColumn(stdData);
stdPred = SLmakeColumn(stdPred);

errorOfmeanPred = SLmakeColumn(errorOfmeanPred);   
errorOfmeanData = SLmakeColumn(errorOfmeanData);   
errorOfstdData = SLmakeColumn(errorOfstdData);   
errorOfstdPred = SLmakeColumn(errorOfstdPred);   


%case meanData and errorData are not input
% if isempty(meanData)==1
%     meanData = nan(size(dataCond,1),1);
% end
% if isempty(errorData)==1
%     errorData = nan(size(dataCond,1),1);
% end
% if isempty(meanPred)==1
%     meanPred = nan(size(dataCond,1),1);
% end
% if isempty(errorPred)==1
%     errorPred = nan(size(dataCond,1),1);
% end


%----------------------------------------
%factors
%----------------------------------------
F.f1.i = dataCond(:,3);
F.f1.nm = 'd';
F.f1.L = unique(F.f1.i);
F.f1.L = sort(F.f1.L,'ascend');
F.f1.n = numel(F.f1.L);

F.f2.i = dataCond(:,2);
F.f2.nm = 'coh';
F.f2.L = unique(F.f2.i);
F.f2.L = sort(F.f2.L,'descend');
F.f2.n = numel(F.f2.L);

F.f3.i = dataCond(:,1);
F.f3.nm = 'Prior std';
F.f3.L = unique(F.f3.i);
F.f3.L = sort(F.f3.L,'descend');
F.f3.n = numel(F.f3.L);

%case bimodal prior
PriorModesUnq = SLuniqpair(priorModes);

%Graphics
F.f2.color = {[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};

F.f2.colorPre = {[0.2 0 0],...
    [0.97 0.2 0],...
    [0.8 0.4 0],...
    [0.3 0.3 0]};

%---------------------------------------------------
%Mean data and predictions
%---------------------------------------------------
scrsz = get(0,'ScreenSize');
fig1 = figure(1);
set(fig1,'color','w','Position',[0 522 480 198]);
ymax = nan(F.f2.n,1);
ymin = nan(F.f2.n,1);
h = nan(F.f2.n,1);
count = 0;

%graphics
marksz = 100;
ftsz = 8;

xCenteredSortedBkp = [];
for j = 1 : F.f2.n
    
    h(j) = subplot(1,F.f2.n,j);
    axis square
    
    for i = 1 : F.f3.n
        
        %data for this conditions
        thisC = F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %recentering to prior
        %--------------------
        %recenter all data relative to prior
        %calculate distance to prior mean 225 deg
        %sort distance
        if isempty(F.f1.i( thisC ))
            sprintf('No data here')
            keyboard
        end
        xCentered = round(SLvectors2signedAngle(F.f1.i( thisC ),225,'polar'));
        
        %express -180? distance to prior as 180? (same but convenient
        %for visualization)
        xCentered(xCentered==-180)=180;
        [xCenteredSorted,IA] = sort(xCentered,'ascend');

        %Predictions
        %(case we want y centered)
        if sum(strcmp(varargin{:},'yCentered'))==1
            
            yCentered = round(SLvectors2signedAngle(meanPred( thisC ),225,'polar'));
            ePredCentered = errorOfmeanPred( thisC );
            
            %express -180? distance to prior as 180? (same but convenient
            %for visualization)
            yCentered(yCentered==-180)=180;
            
        elseif sum(strcmp(varargin{:},'yNotCentered'))==1
            yCentered = meanPred( thisC );
            ePredCentered = errorOfmeanPred( thisC );
        end
        yCenteredSorted = yCentered(IA);
        ePredCenteredSorted = ePredCentered(IA);


        %data
        %(case we want y centered)
        if sum(strcmp(varargin{:},'yCentered'))==1
            
            yDataCentered = round(SLvectors2signedAngle(meanData( thisC ),225,'polar'));
            eDataCentered =  errorOfmeanData( thisC );
            
            %express -180? distance to prior as 180? (same but convenient
            %for visualization)
            yDataCentered(yDataCentered==-180)=180;        
            
        elseif sum(strcmp(varargin{:},'yNotCentered'))==1
            yDataCentered = meanData( thisC );
            eDataCentered = errorOfmeanData( thisC );
        end
        yDataCenteredSorted = yDataCentered(IA);
        eDataCenteredSorted = eDataCentered(IA);
        
        %sort positions and labels for plot such that 225 deg the prior
        %mean is at the center of the plot (when distance = 0).
        %Sorting is done both for x and y data such that they remain
        %aligned across conditions.
        xTickCentered = SLvectors2signedAngle(F.f1.L,225,'polar');
        
        %express -180? distance to prior as 180? (same but convenient
        %for visualization)
        xTickCentered(xTickCentered==-180)=180;
        [~,I] = sort(xTickCentered,'ascend');
        xTickCenteredSorted = xTickCentered(I);
        yTickCenteredSorted = xTickCenteredSorted;
        
        %To plot estimates mean against motion directions on a 
        %linear space (Fig.3A), the raw motion directions and 
        %estimates mean were normalized to vectorial angles from
        %the prior mean (225?) and the x and y axes were centered 
        %at zero (normalized prior mean) by shifting them circularly. 
        %Rotation angles were then labelled according to their raw values
        %on the circle (e.g., 0?, is labelled 225?). Those operations were 
        %performed for both subjects and model mean estimates. The operations 
        %were repeated in figure 3B,C and D but only for the x-axes. A mean 
        %estimate of 33? was calculated for 55? motion direction which is 
        %very far from motion direction on the linear space but actually 
        %close to motion directions on the circular space. We got rid of 
        %this visual artifact by expressing both 55? and 33? as the 
        %counterclockwise distance to prior mean (55? becomes 190? instead 
        %of 170? and 33? becomes 168?). Note that the maximum vectorial angle is  >180?. 
        if j==3 && ~isempty(intersect(xCenteredSorted,180))
            
            %move point at -170? distance to prior at 190? (positive side) 
            %and convert values at x=-170? to positive distance relative to prior
            %to improve visualization
            posNeg170 = xCenteredSorted == -170;
            xCenteredSorted(posNeg170) = 225 - 170 + 360 - 225;
            xCenteredSorted(xCenteredSorted == 180) = - 180;
            yCenteredSorted(posNeg170) = 225 - abs(yCenteredSorted(posNeg170)) + 360 - 225;
            yCenteredSorted(yCenteredSorted == 180) = - 180;
            yDataCenteredSorted(posNeg170) = 225 - abs(yDataCenteredSorted(posNeg170)) + 360 - 225;
            
            %sort
            [xCenteredSorted,I] = sort(xCenteredSorted,'ascend');
            yDataCenteredSorted = yDataCenteredSorted(I);
            yCenteredSorted = yCenteredSorted(I);
            
            %label
            xTickCenteredSorted = xCenteredSorted;
            yTickCenteredSorted = xCenteredSorted; 
        end
        xtickLabel = F.f1.L(I);
        ytickLabel = F.f1.L(I);
        
        
        %error plots
        %-------------
        count = count+1;
%         myPlot1(count) = scatter(xCenteredSorted,yDataCenteredSorted,...
%             marksz,...
%             'MarkerEdgeColor','w',...
%             'MarkerFaceColor',F.f2.color{i},...
%             'displayname',strcat(F.f3.nm,':',num2str(F.f3.L(i))));

%       errorbar for data
        myPlot1(count) = SLerrorbar(xCenteredSorted,yDataCenteredSorted,...
            'yError',eDataCenteredSorted,'Symbol=o',['Color=[' num2str(F.f2.color{i}) ']']);                        
        %erroarea
        %sldrawErrorArea(yDataCenteredSorted,eDataCenteredSorted,xCenteredSorted,'color',F.f2.color{i})
        
%         myPlot2(count) = plot(xCenteredSorted,yCenteredSorted,...
%             'color',F.f2.colorPre{i},...
%             'linewidth',3,...
%             'linestyle','-',...
%             'linesmoothing','on',...
%             'displayName','Bayes');        
%         myPlot2(count) = SLerrorbar(xCenteredSorted,yCenteredSorted,...
%             'yError',ePredCenteredSorted,'Symbol=-',['color=' num2str(F.f2.colorPre{i})]);
        
        %errorarea for model
        sldrawErrorArea(yCenteredSorted,ePredCenteredSorted,xCenteredSorted,'color',F.f2.colorPre{i});
                
        %graphics
        ymax(j,i) = max([yDataCenteredSorted; yCenteredSorted]);
        ymin(j,i) = min([yDataCenteredSorted; yCenteredSorted]);
        
        %backup
        xCenteredSortedBkp = [xCenteredSortedBkp {xCentered}];
    end
    
    %x and ylabel
    if j==1
        ylabel('Mean estimate (deg)','fontsize',ftsz,...
            'FontName','Helvetica',...
            'FontWeight','light',...
            'FontAngle','oblique');
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz,'FontName','Helvetica',...
            'FontWeight','light',...
            'FontAngle','oblique');
    end
    
    [~,pos] = intersect(ytickLabel,[65 145 225 305 25]);

    set(gca,'fontsize',ftsz,'FontName','Helvetica',...
            'FontWeight','light',...
            'FontAngle','oblique',...
            'LineWidth',0.5,...
            'ytick',yTickCenteredSorted(sort(pos)),...
            'yticklabel',ytickLabel(sort(pos)))
    
    %case von Mises prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')
        
        %prior mean marker
        plot([xCenteredSorted(1) xCenteredSorted(end)],[0 0],'b--')
        
        %xtick
        
        [~,pos] = intersect(ytickLabel,[65 145 225 305 25]);

        set(gca,'xtick',xTickCenteredSorted(sort(pos)),...
            'xticklabel',xtickLabel(sort(pos)),...
            'LineWidth',0.5)
                
        %case bimodal prior
        %------------------
    elseif strcmp(priorShape,'bimodalPrior')
        
        for thisidx = 1 : numel(xCenteredSortedBkp)
            numDataPoint(thisidx) = numel(xCenteredSortedBkp{thisidx});
        end
        [~,maxpos] = max(numDataPoint);
        xCenteredSortedFull = xCenteredSortedBkp{maxpos};
        
        %prior mean marker
        thisPriorModes = PriorModesUnq(:);
        for thisidx = 1 : numel(thisPriorModes)
            [~,pos] = intersect(ytickLabel,thisPriorModes(thisidx));
            plot([xCenteredSortedFull(1) xCenteredSortedFull(end)],...
                [yTickCenteredSorted(pos) yTickCenteredSorted(pos)],'b--')
        end
        
        %display prior modes
        [~,pos] = intersect(xtickLabel,thisPriorModes);
        set(gca,'xtick',xTickCenteredSorted(pos),...
            'xticklabel',PriorModesUnq(:),...
            'LineWidth',0.5)
        
        %Warning
    else
        error('(SLdrawModelsPredictionCentered) you need to input prior type')
    end
    
    %selected legend
    %lg = legend(myPlot2(j));
    %legend(lg,'boxoff')
end

%x and ylimits
set(h,'ylim',[min(ymin(:)) max(ymax(:))])
set(h,'xlim',[min(ymin(:)) max(ymax(:))])

%clear up
% SLremoveDeadSpace
% SLConventionUp

%---------------------------------------------------
%std data and predictions
%---------------------------------------------------
fig2 = figure(2);
set(fig2,'color','w','Position',[3 264 480 211])
ymax = nan(F.f2.n,1);
ymin = nan(F.f2.n,1);
h = nan(F.f2.n,1);
count = 0;
for j = 1 : F.f2.n
    
    h(j) = subplot(1,F.f2.n,j);
    axis square
    
    for i = 1 : F.f3.n
        
        %this conditions data
        thisC = F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %transform to angle relative to prior
        %this automatically align prior mean to the center of the plot
        %(distance=0)
        %distance to prior mean
        %sort distance to prior mean
        xCentered = SLvectors2signedAngle(F.f1.i( thisC ),225,'polar');

        %express -180? distance to prior as 180? (same but convenient
        %for visualization)
        xCentered(xCentered==-180)=180;
        
        %center prior mean on x-axis by shift x-axis circularly
        [xCenteredSorted,IA] = sort(xCentered,'ascend');
        
        %move data and pred accordingly
        yDataCentered = stdData( thisC );
        yDataCenteredSorted = yDataCentered(IA);        
        yCentered = stdPred( thisC );
        yCenteredSorted = yCentered(IA);
        
        eDataCentered = errorOfstdData( thisC );
        eDataCenteredSorted = eDataCentered(IA);
        
        ePredCentered = errorOfstdPred( thisC );
        ePredCenteredSorted = ePredCentered(IA);
        
        %labels
        xTickCentered = SLvectors2signedAngle(F.f1.L,225,'polar');
                
        %express -180? distance to prior as 180? to conserve same format as
        %for estimates mean
        xTickCentered(xTickCentered==-180) = 180;
        [xTickCenteredSorted,I] = sort(xTickCentered,'ascend');
        xtickLabel = F.f1.L(I);
        
        %Because of the circular space the y value of the data points far 
        %from prior mean (e.g., -170) look very different from their x value. 
        %But they are actually close to each others. I remove this effect by 
        %converting the y and x values of those points to the distance to the prior 
        %mean on a linear space. Thus y and x are sometimes > or < 180? which 
        %is the max distance to prior mean on a circular space (e.g., -170? move to 190?)
        if j==3 && ~isempty(intersect(xCenteredSorted,180))
            
            %move point at -170? distance to prior at 190? (positive side) 
            %and convert values at x=-170? to positive distance relative to prior
            %to conserve the same format as for estimates mean
            posNeg170 = xCenteredSorted == -170;
            xCenteredSorted(posNeg170) = 225 - 170 + 360 - 225;
            xCenteredSorted(xCenteredSorted == 180) = - 180;
            yCenteredSorted(posNeg170) = 225 - abs(yCenteredSorted(posNeg170)) + 360 - 225;
            yCenteredSorted(yCenteredSorted == 180) = - 180;
            yDataCenteredSorted(posNeg170) = 225 - abs(yDataCenteredSorted(posNeg170)) + 360 - 225;
            
            %sort
            [xCenteredSorted,I] = sort(xCenteredSorted,'ascend');
            yDataCenteredSorted = yDataCenteredSorted(I);
            yCenteredSorted = yCenteredSorted(I);
            
            %label
            xTickCenteredSorted = xCenteredSorted;
            yTickCenteredSorted = xCenteredSorted; 
        end
        
        %plot
        count = count + 1;
        
%         myPlot1(count) = scatter(xCenteredSorted,yDataCenteredSorted,...
%             marksz,...
%             'MarkerEdgeColor','w',...
%             'MarkerFaceColor',F.f2.color{i},...
%             'displayname',strcat(F.f3.nm,':',num2str(F.f3.L(i))));        
%         
%         myPlot2(count) = plot(xCenteredSorted,yCenteredSorted,...
%             'color',F.f2.colorPre{i},...
%             'linewidth',3,...
%             'linestyle','-',...
%             'linesmoothing','on',...
%             'displayName','Bayes');
        
        %errorbar for data
        myPlot1(count) = SLerrorbar(xCenteredSorted,yDataCenteredSorted,...
            'yError',eDataCenteredSorted,'Symbol=o',['Color=[' num2str(F.f2.color{i}) ']']);        
                
        %myPlot2(count) = SLerrorbar(xCenteredSorted,yCenteredSorted,...
        %    'yError',ePredCenteredSorted,'Symbol=-',['color=' num2str(F.f2.colorPre{i})]);
        %errorarea for model
        sldrawErrorArea(yCenteredSorted,ePredCenteredSorted,xCenteredSorted,'color',F.f2.colorPre{i})

        
        %graphics
        xmax(j,i) = max(xCenteredSorted);
        xmin(j,i) = min(xCenteredSorted);
    end
    
    %x and ylimits
    ymin(j,i) = min([yDataCenteredSorted;yCenteredSorted]);
    
    %x and ylabel
    if j==1
        ylabel('Estimate std (deg)','fontsize',ftsz,...
            'FontName','Helvetica',...
            'FontWeight','light',...
            'FontAngle','oblique')
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz,...
            'FontName','Helvetica',...
            'FontWeight','light',...
            'FontAngle','oblique')
    end
    
    set(gca,'fontsize',ftsz,'FontName','Helvetica',...
        'FontWeight','light',...
        'FontAngle','oblique',...
        'LineWidth',0.5,...
        'ytick',0 : 20 : 120,...
        'yticklabel',0 : 20 : 120);
    
        %case von Mises prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')        
        [~,pos] = intersect(xtickLabel,[65 145 225 305 25]);        
        set(gca,'xtick',xTickCenteredSorted(sort(pos)),...
            'xticklabel',xtickLabel(sort(pos)),'LineWidth',0.5,...
            'LineWidth',0.5)        
        %case bimodal prior
        %------------------
    elseif strcmp(priorShape,'bimodalPrior')        
        %display prior modes
        [~,pos] = intersect(xtickLabel,PriorModesUnq(:));
        set(gca,'xtick',xTickCenteredSorted(pos),...
            'xticklabel',PriorModesUnq(:),'LineWidth',0.5)        
        %Warning
    else
        error('(drawMeanPreCentered) you need to input prior type')
    end    
    %selected legend
%     lg = legend(myPlot2(j));
%     legend(lg,'boxoff')
end

%x and ylimits
set(h,'ylim',[0 120])
set(h,'xlim',[min(xmin(:)) - 1  max(xmax(:)) + 1])

