
%author: steeve laquitaine
  %date: 140411
 %usage: SLsortByParamSelectivity(dataTRAINmn,paramValsTrain)
     %e.g., dataTRAINmn is a matrix (m,n)
     %e.g., paramValsTrain is an arraw of parameters (e.g., directions) 
     %with m values.
     
     %dataTRAINmn=rand(67,288);
     %paramValsTrain=randi(360,288,1);
     %SLsortByParamSelectivity(dataTRAINmn,paramValsTrain)

function [dataTRAINmn,VarParamSelectivity,Variable] = SLsortByParamSelectivity(dataTRAINmn,paramValsTrain,VarParamSelectivity,OptionDisp)

%number of rows (e.g., voxels) and columns (trials)
numrow=size(dataTRAINmn,1);

%sort matrix of data by row's preferred parameter value in column
[VarParamSelectivity,Variable]=sort(VarParamSelectivity);
dataTRAINmn=dataTRAINmn(Variable,:);

%draw average data value for single parameter
paramValsTrainUnq=unique(paramValsTrain);
meanMTBOLDoverRep=nan(numrow,numel(paramValsTrainUnq));

%average data over parameter repetitions
for i=1:numel(paramValsTrainUnq);
    meanMTBOLDoverRep(:,i)=mean(dataTRAINmn(:,paramValsTrain==paramValsTrainUnq(i)),2);
end

if strcmp(OptionDisp,'display=on')
    
    %set figure parameters
    screen=get(0,'ScreenSize');
    figure('position',[0.4*screen(3) 0.9*screen(4) 0.4*[screen(3) screen(4)]],'color','w');
    
    %now plot
    for i=1:numel(paramValsTrainUnq);
        plot(1:1:numrow,meanMTBOLDoverRep(:,i),'color',[.5 .5 .5],'linesmooth','on','linewidth',2)
        xlabel('Rows sorted by preferred parameter','fontsize',14)
        ylabel('Data value','fontsize',14)
        xlim([0 numrow])
        ylim([min(meanMTBOLDoverRep(:)) max(meanMTBOLDoverRep(:))])
        set(gca,'fontsize',14,'xtick',1:10:numrow,'xticklabel',VarParamSelectivity(1:10:end))
        box off
        drawnow
        pause(0.7)
    end
    
    %Distribution of voxels preferred directions for MT and LIP. We expect a
    %uniform representation of preferred directions across voxels.
    figure('color','w','position',[0.4*screen(3) 0.05*screen(4) 0.4*[screen(3) screen(4)]], 'color','w')
    hist(VarParamSelectivity)
    set(get(gca,'child'),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
    box off
    title('MT voxels')
    xlabel('Preferred parameter')
    ylabel('Number of rows')
    box off
    
    %draw raw data sorted by parameter selectivity (column)
    figure('color','w');
    imagesc(dataTRAINmn);
    xlabel('Directions (degrees)','fontsize',14);
    ylabel('Voxels by direction selectivity (degrees)')
    set(gca,'ytick',1:10:numel(VarParamSelectivity),'yticklabel',VarParamSelectivity(1:10:end))
    set(gca,'xtick',1:10:numel(paramValsTrain),'xticklabel',paramValsTrain(1:10:end))
    fprintf('%12s \n','done')
end

