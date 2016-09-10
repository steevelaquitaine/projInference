

%slgetprevDirEffect.m
%
% author: steeve laquitaine
%purpose: get all the trials in which the currently displayed direction is 
%         in between the previous direction and the prior mean and plot the
%         average estimates to test whether they are biased toward the
%         previous direction or the prior mean.
%
%  usage: 
%
%       [meanDispDirInfo,meanEstimateInfo,meanPrevDirInfo,FA1,FA2] = slgetprevDirEffect(databank,FA,varargin)
%
%note : databank can be for one subject or many
%
%
%Description: 
%
%   the data could be explained by two mechanisms (long or short-term memory)
%
%   Prior effect (trial history)
%   ----------------------------
%   Subjects learnt the mean of the experimental prior and
%   use it to guess the noisy motion direction in the current trial.
%   Bayesian and Competition mechanisms are both consistent with this
%   mechanism as those mechanism represent and use the prior
%   (mean and strength).
%
%   Recency effect (just previous trial)
%   ------------------------------------
%   Subjects are simply biased by the past trial's motion
%   direction and do not learn the prior.
%
%So, are data attracted toward the mean the experimental prior (prior
%effect) or toward the past trial's motion direction (recency effect)?


function [meanDispDirInfo,meanEstimateInfo,meanPrevDirInfo,FA1,FA2] = slgetprevDirEffect(databank,FA,varargin)

fprintf('%s \n','(analyses) Now testing for "sequential effect"...')

%(Case von Mises prior experiment)
%--------------------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')           
        
        %Get the coordinates of the estimated directions (actual code)
        dir.estiCoor = cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
        
        %Get the coordinates of the displayed directions
        r = 1;
        dir.dispCoor = polar2cartesian(FA.g3.thisT,r);
        
        %Get the coordinates of the prior mean
        dir.priormeanCoor = polar2cartesian(225,r);
        
        %Get the directions that were displayed between the mean of the prior and
        %the past trial's directions.
        %e.g., 225? (prior mean) < 180? (current) < 135? (previous)
        %Calculate the average direction displayed on current trials, on past
        %trial and look at where average estimate in current trials lean
        %towards: prior or past trial?
        %Look at each factor 1
        meanDispDir = num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        meanPrevDir = num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        meanEstimate = num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        
        %calculate and draw
        for j = 1 : FA.g1.lvlsnb
            
            %find this factors' trials
            trialsThisFA1 = FA.g1.thisT==FA.g1.lvlsnm(j);
            
            %look at factor 2
            for k = 1 : FA.g2.lvlsnb
                
                %find this factor's trials
                trialsThisFA2 = FA.g2.thisT==FA.g2.lvlsnm(k);
                
                %get positions of current trials with this condition
                posCurrentThisCondition=find(trialsThisFA1 & trialsThisFA2);
                
                %get positions of the trials that preceded
                posPrevious=posCurrentThisCondition-1;
                
                %remove first trials because there is not preceding trial
                posPrevious(posCurrentThisCondition==1)=[];
                posCurrentThisCondition(posCurrentThisCondition==1)=[];
                
                %get current and previous direction
                dir.dispCoorCurThisCondition=dir.dispCoor(posCurrentThisCondition,:);
                dir.dispCoorPrev=dir.dispCoor(posPrevious,:);
                
                %get the current directions displayed in our condition
                for i = 1 : size(dir.dispCoorCurThisCondition,1)
                    
                    %When the current and the previous directions were displayed
                    %counterclockwise to prior's mean (i.e.,0< angle(prior,current)<180)
                    angle.PriortoCur(i) = SLvectors2signedAngle(dir.priormeanCoor,...
                        dir.dispCoorCurThisCondition(i,:),'cartesian');
                    
                    %angle between current and previous direction
                    angle.CurtoPrev(i) = SLvectors2signedAngle(dir.dispCoorCurThisCondition(i,:),...
                        dir.dispCoorPrev(i,:),'cartesian');
                    
                    %get current trials when the sequence is: past < current < prior
                    %case current counterclockwise to prior mean
                    %case previous counterclockwise to current for this condition
                    if  0<angle.PriortoCur(i) && angle.PriortoCur(i)<180 && 0<angle.CurtoPrev(i) && angle.CurtoPrev(i)<180
                        
                        %store the position of those trials
                        trialCounterClock{j,k}(i) = posCurrentThisCondition(i);
                    else
                        trialCounterClock{j,k}(i) = NaN;
                    end
                end
                
                %now get rid of NaN
                trialCounterClock{j,k}(isnan(trialCounterClock{j,k}))=[];
                
                
                %-------------------------------------------------
                %averages current, previous directions & estimates
                %-------------------------------------------------
                %displayed direction in current trial
                meanDispDirInfo{j,k} = SLvectorStat(dir.dispCoor(trialCounterClock{j,k},:));
                meanDispDir{j,k} = meanDispDirInfo{j,k}.coord.mean;
                
                %estimated direction in current trial
                meanEstimateInfo{j,k} = SLvectorStat(dir.estiCoor(trialCounterClock{j,k},:));
                meanEstimate{j,k} = meanEstimateInfo{j,k}.coord.mean;
                
                %displayed direction in previous trial
                meanPrevDirInfo{j,k} = SLvectorStat(dir.dispCoor(trialCounterClock{j,k}-1,:));
                meanPrevDir{j,k} = meanPrevDirInfo{j,k}.coord.mean;
            end
        end

        %save associated factors
        FA1 = repmat(FA.g1.lvlsnm,1,FA.g2.lvlsnb);
        FA2 = repmat(FA.g2.lvlsnm',FA.g1.lvlsnb,1);        
        
        %draw all conditions
        meanPriorCoord=polar2cartesian(225,10);
        count=0;
        
        %graphics
        fig=figure('color','w');
        ax1=[];
        FA1nmlabel=[];
        FA2nmlabel=[];
        
        %look at each condition
        for k=1:FA.g2.lvlsnb
            
            for j=1:FA.g1.lvlsnb
                
                %draw Estimates density
                %----------------------
                %get all estimates of current trials directions in this condition
                CurEstimateThisCondition=meanEstimateInfo{j,k}.deg.all;
                
                %get number of estimates
                numCurEstimateThisCondition=length(CurEstimateThisCondition);
                
                %count occurrence of each estimated direction
                x = 5 : 10 : 355;
                probaCurEstimateThisCondition = histc(CurEstimateThisCondition,x)/numCurEstimateThisCondition;
                
                %initialize axes and axes' legends
                count = count + 1;
                ax1=[ax1 subplot(numel(FA.g2.lvlsnm),numel(FA.g1.lvlsnm),count)];
                FA1nmlabel=[FA1nmlabel FA.g1.lvlsnm(j)];
                FA2nmlabel=[FA2nmlabel FA.g2.lvlsnm(k)];
                
                %get the maximum of the estimate density to scale mean estimate
                %vector
                maxData=max(probaCurEstimateThisCondition);
                
                %increase the polar radius to get a bit of space for visibility
                radius=1.1*maxData;
                %SLpolar(de2r(x,[]),probaCurEstimateThisCondition',[0.7 0 0],'xSpokes',...
                %    'SpokesTickSteps',8,'area',...
                %    'facecolor',[1 0.5 0.5],...
                %    'edgecolor',[1 0.5 0.5],...
                %    'linestyle','-')
                
                SLpolar(SLde2r(x,[]),probaCurEstimateThisCondition',[0.7 0 0],...
                    'edgecolor',[1 0.5 0.5],...
                    'linestyle','-')                
                
                %draw the mean estimate vector
                %-----------------------------
                %get lenth of vector
                Vectnorm(j,k) = SLcalculateVectNorm(meanEstimate{j,k}(1),meanEstimate{j,k}(2));
                
                %this scales the arrow to the size of the radius
                putOncircle=Vectnorm(j,k)/maxData;
                
                %draw mean estimated direction
                hold on
                arrow([0 0],[meanEstimate{j,k}(1) meanEstimate{j,k}(2)]/putOncircle,...
                    'Length',2,'BaseAngle',[],'Width',1,...
                    'facecolor',[0.5 0 0],...
                    'edgecolor',[0.5 0 0],...
                    'linestyle','-')
                
                %draw current direction
                %----------------------
                %indicate the mean of the current direction
                %get all current trials displayed directions in this condition
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanDispDir{j,k}(1),meanDispDir{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanDispDir{j,k}(1)/putOncircle,meanDispDir{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[.2 .2 .2])
                
                %draw previous direction
                %-----------------------
                %indicate the mean of the current direction
                %get all current trials displayed directions in this condition
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanPrevDir{j,k}(1),meanPrevDir{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanPrevDir{j,k}(1)/putOncircle,meanPrevDir{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[.7 .7 .7])
                
                
                %draw prior mean
                %---------------
                %indicate the mean of the prior
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanPriorCoord(1),meanPriorCoord(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanPriorCoord(1)/putOncircle,meanPriorCoord(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[0 0.7 1])
            end
        end
        
        %legend the axes
        for i=1:numel(ax1)
            title(ax1(i),[num2str(100*FA1nmlabel(i)),'% ',FA.g1.nm,' - ',...
                num2str(FA2nmlabel(i)),' degrees ',FA.g2.nm],'fontsize',8,...
                'fontweight','Bold')
        end
        
        %remove dead spaces
        SLremoveDeadSpace
        
    end
end

% 
% %(Case bimodal prior experiment)
% %-------------------------------
% if sum(strcmp(varargin{1},'experiment'))
%     if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')%
%         
%         %Get the coordinates of the estimated directions (actual code)
%         dir.estiCoor=cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
%         
%         %Test the code (all tests worked)
%         %----------------------------------------------------------------------
%         %%simulate when subjects estimate the current direction: arrow should point
%         %current direction if code correct
%         %dir.estiCoor=polar2cartesian(FA.g3.thisT,1);
%         
%         %%simulate when subjects estimate the previous direction: arrow should point
%         %previous direction if code correct
%         %dir.estiCoor=polar2cartesian([NaN;FA.g3.thisT(1:end-1)],1);
%         %----------------------------------------------------------------------
%         
%         %set the radius of the circle but the true size is 2.5 degrees
%         %(just for visualization purpose)
%         r=1;
%         
%         %Get the coordinates of the displayed directions
%         dir.dispCoor=polar2cartesian(FA.g3.thisT,r);
%         
%         %counterclock mode at each trial
%         %when the distance between the first and second mode is <0 the
%         %first mode is more counterclockwise than the second and when it is
%         %>0, the second mode is more clockwise than the first. Now get the
%         %coordinates of the counterclock mode (left)
%         distModes=SLvectors2signedAngle(FA.g2.thisT(:,1),FA.g2.thisT(:,2));
%         cclockMode(distModes<0)=FA.g2.thisT(distModes<0,1);
%         cclockMode(distModes>0)=FA.g2.thisT(distModes>0,2);
%         dir.priorcclockModeCoor=polar2cartesian(cclockMode,r);
%         
%         %Now get the coordinates of the clock mode (right)
%         clockMode(distModes<0)=FA.g2.thisT(distModes<0,2);
%         clockMode(distModes>0)=FA.g2.thisT(distModes>0,1);
%         dir.priorclockModeCoor=polar2cartesian(clockMode,r);
%         
%         %Get the directions that were displayed between the mean of the prior and
%         %the past trial's directions.
%         %e.g., 225? (prior mean) < 180? (current) < 135? (previous)
%         %Calculate teh average direction displayed on current trials, on past
%         %trial and look at where average estimate in current trials lean
%         %towards: prior or past trial?
%         %Look at each factor 1
%         meanDispDir=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
%         meanPrevDir=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
%         meanEstimate=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
%         
%         %case factor 1 is the prior modes
%         if strcmp(FA.g1.nm,'priormodes')
%             FA.g1.thisTSeqAna=FA.g1.distance_mode;
%             %case not prior mode
%         else
%             FA.g1.thisTSeqAna=FA.g1.thisT;
%         end
%         
%         %case factor 2 is the prior modes
%         if strcmp(FA.g2.nm,'priormodes')
%             FA.g2.thisTSeqAna=FA.g2.distance_mode;
%             %case not prior mode
%         else
%             FA.g2.thisTSeqAna=FA.g2.thisT;
%         end
%         
%         %case factor 3 is the prior modes
%         if strcmp(FA.g3.nm,'priormodes')
%             FA.g3.thisTSeqAna=FA.g3.distance_mode;
%             %case not prior mode
%         else
%             FA.g3.thisTSeqAna=FA.g3.thisT;
%         end
%         
%         %calculate and draw
%         for j=1:FA.g1.lvlsnb
%             
%             %find this factors' trials
%             trialsThisFA1=FA.g1.thisTSeqAna==FA.g1.lvlsnm(j);
%             
%             %look at factor 2
%             for k=1:FA.g2.lvlsnb
%                 
%                 %find this factor's trials
%                 trialsThisFA2=FA.g2.thisTSeqAna==FA.g2.lvlsnm(k);
%                 
%                 %get positions of current trials with this condition
%                 posCurrentThisCondition=find(trialsThisFA1 & trialsThisFA2);
%                 
%                 %get positions of the trials that preceded
%                 posPrevious=posCurrentThisCondition-1;
%                 
%                 %remove first trials because there is not a preceding trial
%                 posPrevious(posCurrentThisCondition==1)=[];
%                 posCurrentThisCondition(posCurrentThisCondition==1)=[];
%                 
%                 %get current trial's motion direction
%                 dir.dispCoorCurThisCondition=dir.dispCoor(posCurrentThisCondition,:);
%                 
%                 %get counterclock mode in the current trial
%                 dir.priorcclockModeCoorThisCondition=dir.priorcclockModeCoor(posCurrentThisCondition,:);
%                 
%                 %get the previous motion direction
%                 dir.dispCoorPrev=dir.dispCoor(posPrevious,:);
%                 
%                 
%                 
%                 %get the current directions
%                 for i=1:size(dir.dispCoorCurThisCondition,1)
%                     
%                     %When current and previous directions were displayed
%                     %counterclockwise to the most counterclockwise mode of
%                     %the bimodal prior (i.e.,0< angle(left prior mode,current)<180)
%                     angle.PriortoCur(i)=SLvectors2signedAngle(dir.priorcclockModeCoorThisCondition(i,:),...
%                         dir.dispCoorCurThisCondition(i,:));
%                     
%                     %angle between current and previous direction
%                     angle.CurtoPrev(i)=SLvectors2signedAngle(dir.dispCoorCurThisCondition(i,:),...
%                         dir.dispCoorPrev(i,:));
%                     
%                     %get current trials when the sequence is: past < current < prior
%                     %case current counterclockwise to prior mean
%                     %case previous counterclockwise to current for this condition
%                     if  0<angle.PriortoCur(i) && angle.PriortoCur(i)<180 && 0<angle.CurtoPrev(i) && angle.CurtoPrev(i)<180
%                         
%                         %store the position of those trials
%                         trialCounterClock{j,k}(i)=posCurrentThisCondition(i);
%                     else
%                         trialCounterClock{j,k}(i)=NaN;
%                     end
%                 end
%                 
%                 %now get rid of NaN
%                 trialCounterClock{j,k}(isnan(trialCounterClock{j,k}))=[];
%                 
%                 %averages current, previous directions & estimates
%                 %-------------------------------------------------
%                 %displayed direction in current trial
%                 meanDispDirInfo{j,k}=SLvectorStat(dir.dispCoor(trialCounterClock{j,k},:));
%                 meanDispDir{j,k}=meanDispDirInfo{j,k}.coord.mean;
%                 
%                 %estimated direction in current trial
%                 meanEstimateInfo{j,k}=SLvectorStat(dir.estiCoor(trialCounterClock{j,k},:));
%                 meanEstimate{j,k}=meanEstimateInfo{j,k}.coord.mean;
%                 
%                 %counterclock prior mode in current trial
%                 meancclockPmodeInfo{j,k}=SLvectorStat(dir.priorcclockModeCoor(trialCounterClock{j,k},:));
%                 meancclockPmode{j,k}=meancclockPmodeInfo{j,k}.coord.mean;
%                 
%                 %clock prior mode in current trial
%                 meanclockPmodeInfo{j,k}=SLvectorStat(dir.priorclockModeCoor(trialCounterClock{j,k},:));
%                 meanclockPmode{j,k}=meanclockPmodeInfo{j,k}.coord.mean;
%                 
%                 %displayed direction in previous trial
%                 meanPrevDirInfo{j,k}=SLvectorStat(dir.dispCoor(trialCounterClock{j,k}-1,:));
%                 meanPrevDir{j,k}=meanPrevDirInfo{j,k}.coord.mean;
%             end
%         end
%         
%         %draw all conditions
%         count=0;
%         
%         %graphics
%         figure('color','w');
%         ax1=[];
%         FA1nmlabel=[];
%         FA2nmlabel=[];
%         
%         %look at each condition
%         for k=1:FA.g2.lvlsnb
%             
%             for j=1:FA.g1.lvlsnb
%                 
%                 %draw Estimates density
%                 %----------------------
%                 %get all estimates of current trials directions in this condition
%                 CurEstimateThisCondition=meanEstimateInfo{j,k}.deg.all;
%                 
%                 %get number of estimates
%                 numCurEstimateThisCondition=length(CurEstimateThisCondition);
%                 
%                 %count occurrence of each estimated direction
%                 x=5:10:355;
%                 probaCurEstimateThisCondition=histc(CurEstimateThisCondition,x)/numCurEstimateThisCondition;
%                 
%                 %initialize axes and axes' legends
%                 count=count+1;
%                 ax1=[ax1 subplot(numel(FA.g2.lvlsnm),numel(FA.g1.lvlsnm),count)];
%                 FA1nmlabel=[FA1nmlabel FA.g1.lvlsnm(j)];
%                 FA2nmlabel=[FA2nmlabel FA.g2.lvlsnm(k)];
%                 
%                 %get the maximum of the estimate density to scale mean estimate
%                 %vector
%                 maxData=max(probaCurEstimateThisCondition);
%                 
%                 %increase the polar radius to get a bit of space for visibility
%                 radius=1.1*maxData;
%                 
%                 %if we found data when motion direction is in between counterclock
%                 %prior mode and previous trial
%                 if ~isnan(CurEstimateThisCondition)
%                     
%                     %plot polar
%                     SLpolar(de2r(x,[]),probaCurEstimateThisCondition,[0.7 0 0],'xSpokes',...
%                         'SpokesTickSteps',4,'area',...
%                         'facecolor',[1 0.5 0.5],...
%                         'edgecolor',[1 0.5 0.5])
%                     
%                     %draw the mean estimate vector
%                     %get lenth of vector
%                     Vectnorm(j,k)=SLcalculateVectNorm(meanEstimate{j,k}(1),meanEstimate{j,k}(2));
%                     
%                     %this scales the arrow to the size of the radius
%                     putOncircle=Vectnorm(j,k)/maxData;
%                     
%                 else
%                     %plot blank polar and null variable
%                     SLpolar(de2r(x,[]),0.1*ones(numel(de2r(x,[])),1),[1 1 1],'xSpokes',...
%                         'SpokesTickSteps',4,'area',...
%                         'facecolor',[1 1 1],...
%                         'edgecolor',[1 1 1])
%                     meanEstimate{j,k}=[0 0];
%                     putOncircle=1;
%                 end
%                 
%                 %draw mean estimated direction
%                 hold on
%                 arrow([0 0],[meanEstimate{j,k}(1) meanEstimate{j,k}(2)]/putOncircle,...
%                     'Length',2,'BaseAngle',[],'Width',1,...
%                     'facecolor',[0.5 0 0],...
%                     'edgecolor',[0.5 0 0])
%                 
%                 %draw current direction
%                 %----------------------
%                 %indicate the mean of the current direction
%                 %get all current trials displayed directions in this condition
%                 %get lenth of vector
%                 Vectnorm(j,k)=SLcalculateVectNorm(meanDispDir{j,k}(1),meanDispDir{j,k}(2));
%                 
%                 %scale its length to the size of the radius
%                 putOncircle=Vectnorm(j,k)/radius;
%                 hold on
%                 plot(meanDispDir{j,k}(1)/putOncircle,meanDispDir{j,k}(2)/putOncircle,'.',...
%                     'MarkerSize',20,...
%                     'color',[.2 .2 .2])
%                 
%                 %draw previous direction
%                 %-----------------------
%                 %indicate the mean of the current direction
%                 %get all current trials displayed directions in this condition
%                 %get lenth of vector
%                 Vectnorm(j,k)=SLcalculateVectNorm(meanPrevDir{j,k}(1),meanPrevDir{j,k}(2));
%                 
%                 %scale its length to the size of the radius
%                 putOncircle=Vectnorm(j,k)/radius;
%                 hold on
%                 plot(meanPrevDir{j,k}(1)/putOncircle,meanPrevDir{j,k}(2)/putOncircle,'.',...
%                     'MarkerSize',20,...
%                     'color',[.7 .7 .7])
%                 
%                 
%                 %draw average prior modes
%                 %------------------------
%                 %indicate the counterclock mode of the prior
%                 %get lenth of vector
%                 Vectnormcclock(j,k)=SLcalculateVectNorm(meancclockPmode{j,k}(1),meancclockPmode{j,k}(2));
%                 
%                 %scale its length to the size of the radius
%                 putOncircle=Vectnormcclock(j,k)/radius;
%                 hold on
%                 plot(meancclockPmode{j,k}(1)/putOncircle,meancclockPmode{j,k}(2)/putOncircle,'.',...
%                     'MarkerSize',20,...
%                     'color',[0 0.7 1])
%                 
%                 %indicate the clock mode of the prior
%                 Vectnormclock(j,k)=SLcalculateVectNorm(meanclockPmode{j,k}(1),meanclockPmode{j,k}(2));
%                 putOncircle=Vectnormclock(j,k)/radius;
%                 hold on
%                 plot(meanclockPmode{j,k}(1)/putOncircle,meanclockPmode{j,k}(2)/putOncircle,'.',...
%                     'MarkerSize',20,...
%                     'color',[0 0.7 1])
%             end
%         end
%         
%         %legend the axes
%         for i=1:numel(ax1)
%             title(ax1(i),[num2str(100*FA1nmlabel(i)),'% ',FA.g1.nm,' - ',...
%                 num2str(FA2nmlabel(i)),' degrees (mode1-mode2) ',FA.g2.nm],'fontsize',8,...
%                 'fontweight','Bold')
%         end
%         
%         %remove dead spaces
%         SLremoveDeadSpace(0.05)
%     end
% end
