
%slplotMeanAndVarBehavior.m
%
%
% author: steeve laquitaine
%   date: 160425
%purpose: plot mean and variability by coherence and prior
%         conditions
%  usage:
%
%e.g., 
%     %database  
%     d = SLMakedatabank({'sub01'},'dataPath','~/data/dataPsychophy/proj01_priorStrength/','experiment','vonMisesPrior');
%     %sort data by variables
%     [fa,estimates] = slInitFAs(d,{'StimStrength','Pstd','FeatureSample'})
%     %visualize
%     [figs,GenStat] = slplotMeanAndVarBehavior(d,fa,estimates,'experiment','vonMisesPrior')

function [figs,GenStat] = slplotMeanAndVarBehavior(databank,FA,est_dir,varargin)

%init figure
fig1.hdle = figure('color','w');
fig1.nm = ['myAnalysis','_mean'];
set(fig1.hdle, 'units', 'centimeters', 'pos', [0 0 30 10])

%graphic parameters
FA.g2.color = {[0.5 0 0],...
    [1 0 0],...
    [1 0.5 0],...
    [0.65 0.65 0]};

%for fit
FA.g2.colorfit = {[0.7 0.2 0],...
    [1 0.2 0],...
    [1 0.7 0],...
    [0.7 0.7 0]};

%(case bimodal prior) get unique prior modes
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
        priormodes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
        uniquePriormodes = unique(sort(priormodes,2), 'rows');
    end
end

%behavior average (mean)
%-----------------------
fprintf('%s \n','(slplotMeanAndVarBehavior) Now plotting the mean of the estimated data...')
for k = 1:FA.g1.lvlsnb
    subplot(1,FA.g1.lvlsnb,k)
    axis square
    hold all
    
    %(case Gaussian prior) indicate the priors' mean
    if sum(strcmp(varargin{1},'experiment'))
        if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
            plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],...
                [databank.data{1,7} databank.data{1,7}],...
                'b:',...
                'linewidth',2);
        end
    end
    
    %(case bimodal prior) indicate the priors' modes
    if sum(strcmp(varargin{1},'experiment'))
        if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
            
            %plot the modes of each prior
            for z=1:size(uniquePriormodes,1)
                
                %get this prior's two modes and plot them
                thesemodes=repmat(uniquePriormodes(z,:),2,1);
                plot(repmat([FA.g3.lvlsnm(1),FA.g3.lvlsnm(end)]',1,2),thesemodes,...
                    ':',...
                    'color',FA.g2.color{z},...
                    'linewidth',1);
            end
        end
    end
    
    %indicate ideal estimations
    plot(FA.g3.lvlsnm,FA.g3.lvlsnm','k:','linewidth',2);
    
    %plot estimates stats for each experimental condition
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            
            %data
            fig1.datamean.deg(j,i,k)=est_dir{j,i,k}.deg.mean;
            fig1.datasem(j,i,k)=est_dir{j,i,k}.deg.sem;
            fig1.datamean.coord{j,i,k}=est_dir{j,i,k}.coord.mean;
            
            %factors g1,g2,g3
            fig1.dataInfog1{j,i,k}=FA.g1.lvlsnm(k);
            fig1.dataInfog2{j,i,k}=FA.g2.lvlsnm(j);
            fig1.dataInfog3{j,i,k}=FA.g3.lvlsnm(i);
            
            %motion directions in cartesian coordinates
            r=2.5;
            %fig1.displayed.coord{j,i,k}=polar2cartesian(fig1.dataInfog3{j,i,k},r);
            
            %Geometric method
            %Calculate distance between displayed and estimated dirs.
            %fig1.Disp_dataDistance{j,i,k}=SLvectors2signedAngle(fig1.datamean.coord{j,i,k}, fig1.displayed.coord{j,i,k} );
            
            %Select raw data
            fig1.datamean.deg(j,i,k)=fig1.datamean.deg(j,i,k);
        end
        
        %Linear fit
        sllinefit(FA.g3.lvlsnm',fig1.datamean.deg(j,:,k),1:1:360,FA.g2.colorfit{j});
        
        %plot mean data
        scatter(FA.g3.lvlsnm,fig1.datamean.deg(j,:,k),...
            100,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',FA.g2.color{j},...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))))
        
        %plot sem
        SLerrorbar(FA.g3.lvlsnm,fig1.datamean.deg(j,:,k),...
            'yError',fig1.datasem(j,:,k),...
            'MarkerSize=1',['Color=',num2str(FA.g2.color{j})])
    end
    
    %graphics
    %--------
    set(gca,...
        'xtick', 0:45:360,'xticklabel',0:45:360,...
        'ytick', 0:45:360,'yticklabel',0:45:360,...
        'fontsize',13);
    xlim([0 360]);
    ylim([0 360]);
    
    %set labels
    if k==1
        ylabel('Estimated direction (deg)','fontsize',14);
    end
    xlabel('Displayed direction (deg)','fontsize',14);
    
    %remove dead space
    drawPublishAxis
end
SLremoveDeadSpace



%behavior variability (std)
%-------------------------
%Draw
fprintf('%s \n','(slplotMeanAndVarBehavior) Now plotting the std of the data...')
fig2.hdle=figure('color','w');
fig2.nm='MyAnalysis_std';
set(fig2.hdle,'units','centimeters','pos',[0 0 30 10])

%Set the number of bootstrap
numboot=50;
tic
maxPlot=[];
for k=1:FA.g1.lvlsnb
    subplot(1,3,k)
    axis square
    hold all
    
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            
            %data
            fig2.datastd(j,i,k)=est_dir{j,i,k}.deg.std;
            fig2.datanb(j,i,k)=est_dir{j,i,k}.num;
            
            %bootstrap data (density of std)
            %initialize the data to bootstrap
            data2boot=[];
            fig2.databootstdSeries{j,i,k}=[];
            
            %collect the data to bootstrap (e.g. data at each displayed dir)
            fig2.datacoord(j,i,k)=est_dir{j,i,k}.coord;
            data2boot.data=fig2.datacoord(j,i,k).all;
            
            %create samples
            for jj=1:numboot
                
                %count the data
                data2boot.nb=size(data2boot.data,1);
                
                %if data are observed in x and number of sample >=2
                if isempty(data2boot.data)==0&&data2boot.nb>=2
                    
                    %bootstrap the data in x
                    %sample randomly a data value from x 'data2boot.nb' times
                    for iii=1:data2boot.nb
                        jit=ceil(rand*data2boot.nb);
                        data_boottmp{j,i,k}{jj}(iii,:)=data2boot.data(jit,:);
                        
                        %store for fig2
                        fig2.databoot{j,i,k}{jj}(iii,:)=data_boottmp{j,i,k}{jj}(iii,:) ;
                        fig2.databootnb{j,i,k}{jj}(iii,:)=data2boot.nb;
                    end
                    
                    %calculate the std for each bootstrap
                    fig2.databootstd{j,i,k}{jj}=SLstatcircular(fig2.databoot{j,i,k}{jj});
                    
                    %collect the std in a vector
                    fig2.databootstdSeries{j,i,k}=[fig2.databootstdSeries{j,i,k} fig2.databootstd{j,i,k}{jj}.deg.std];
                    
                    %if there are no data in x
                else
                    %fill in a nan
                    data_boottmp{j,i,k}{jj}=nan;
                    
                    %store for fig2
                    fig2.databoot{j,i,k}{jj}=data_boottmp{j,i,k}{jj} ;
                    
                    %there are no data
                    fig2.databootnb{j,i,k}{jj}=0 ;
                    
                    %nan the std
                    fig2.databootstd{j,i,k}{jj}=nan;
                    
                    %nan the vector of std
                    fig2.databootstdSeries{j,i,k}=nan;
                end
            end
            
            %calculate the std over bootstrapped samples
            fig2.databootMeanOfstd{j,i,k}=nanmean(fig2.databootstdSeries{j,i,k});
            fig2.databootStdOfstd{j,i,k}=std(fig2.databootstdSeries{j,i,k});
            
            %factors g1,g2,g3
            fig2.dataInfog1{j,i,k}=FA.g1.lvlsnm(k);
            fig2.dataInfog2{j,i,k}=FA.g2.lvlsnm(j);
            fig2.dataInfog3{j,i,k}=FA.g3.lvlsnm(i);
        end
        
        %draw variability
        myerrorbar(FA.g3.lvlsnm,[fig2.databootMeanOfstd{j,:,k}], 'yError',...
            [fig2.databootStdOfstd{j,:,k}],...
            'Symbol=.',...
            'LineWidth=1',...
            ['Color=',num2str(FA.g2.color{j})],...
            ['MarkerFacecolor=', num2str(FA.g2.color{j})],...
            'MarkerSize=1',...
            'MarkerEdgeColor=w');
        
        scatter(FA.g3.lvlsnm,[fig2.databootMeanOfstd{j,:,k}],100,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',FA.g2.color{j},...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))))
        
        h(k)=gca;
        maxPlot=[maxPlot max(sum([[fig2.databootMeanOfstd{j,:,k}]; [fig2.databootStdOfstd{j,:,k}]]))];
    end
    
    %graphics
    %--------
    %x-axis
    set(gca,'xtick',0:45:360,'xticklabel',0:45:360,...
        'fontsize',13)
    set(gca,'ytick',0:10:100,'yticklabel',0:10:100,'fontsize',13);
    xlim([0 360]);
    ylim([0 100]);
    
    %set labels of axes
    if k==1
        ylabel('Std of estimated direction (degrees)','fontsize',11);
    end
    xlabel({'Displayed direction (degrees)','(degrees)'},'fontsize',11);
    
    %make clean
    drawPublishAxis
end
SLremoveDeadSpace
% for k=1:FA.g1.lvlsnb
%     ylim(h(k),[0 max(maxPlot)])
%     set(h(k),'ytick',1:10:max(maxPlot),'yticklabel',1:10:max(maxPlot))
% end
%Draw raw data
% fprintf('%s \n','(drawRawEstDir) Now plotting individual estimates in polar...')
% fig4=drawRawEstDir(fig2);

%backup figures informations
figs = {fig1, fig2};

%case Ancova
if strcmp(varargin,'Ancova')
            
    %%Stat of ANCOVA: effect of FA1(e.g., coherence) on the slope of data=f(FA3)).
    %%----------------------------------------------------------------------
    %%Set x: (e.g.,FA3 - displayed dirs)
    %x=[];
    %x=reshape(repmat(FA.g3.lvlsnm , FA.g2.lvlsnb*FA.g1.lvlsnb,1), [],1);
    %
    %%Set y: (data; e.g., estimated dirs)
    %%Loop over the levels of Factor 1
    %y=[];
    %for k=1:FA.g1.lvlsnb
    %    ytmp=[];
    %    ytmp=reshape(fig1.datamean.deg(:,:,k)',[],1);
    %    y=[y; ytmp];
    %end
    %
    %%Set g the groups (e.g., FA1 - coherences)
    %g=[];
    %g=reshape(repmat(1: FA.g1.lvlsnb, FA.g3.lvlsnb*FA.g2.lvlsnb,1), [],1);
    %
    %%Run the Ancova
    %Ancova_effect_of_FA1=statAncova(x, y, g);
    %display('Your Ancova is being performed...')
    %
    %%%Report the statistics of the Ancova
    %%%Degrees of freedom
    %%%Within groups  (e.g.,=f(nb of subjects))
    %%dfwgT1  =GenStat{2,1}{2,2}{5,2};
    %%%Between groups  (e.g.,=f(nb of group))
    %%dfbgT1  =GenStat{2,1}{2,2}{4,2};
    %%%F-value for the interaction effect
    %%FvalueT1=GenStat{2,1}{2,2}{4,5};
    %%%p-value of the tested interaction effect
    %%pvalueT1=GenStat{2,1}{2,2}{4,6};
    %%GenStat2ReportT1={dfbgT1, dfwgT1, FvalueT1, pvalueT1};
    %%display(GenStat2ReportT1)
    
    %
    %----------------------------------------------------------------------
    %Stat of ANCOVA: effect of FA1(e.g., prior's width) on the slope of data=f(FA3)).
    %----------------------------------------------------------------------
    %Set x: (e.g.,FA3 - displayed dirs)
    x=[];
    x=reshape(repmat(FA.g3.lvlsnm , FA.g2.lvlsnb*FA.g1.lvlsnb,1), [],1);
    
    %Set y: (data; e.g., estimated dirs)
    %Loop over the levels of Factor 1
    y=[];
    for k=1:FA.g1.lvlsnb
        ytmp=[];
        ytmp=reshape(fig1.datamean.deg(:,:,k)',[],1);
        y=[y; ytmp];
    end
    
    %Set g the groups (e.g., FA1 - coherences)
    g=[];
    g=reshape(repmat(1: FA.g2.lvlsnb, FA.g3.lvlsnb*FA.g1.lvlsnb,1), [],1);
    
    %check enough groups
    if length(unique(g))<2
        GenStat = [];
        sprintf('(slplotMeanAndVarBehavior) Not enough groups for ANCOVA')
        return
    end
    
    %Run the Ancova
    Ancova_effect_of_FA2 = slStatAncova(x, y, g);
    fprintf('%s \n','(slplotMeanAndVarBehavior) running Ancova...')
    
    %%Report the statistics of the Ancova
    %%Degrees of freedom
    %%Within groups  (e.g.,=f(nb of subjects))
    %dfwgT1  =GenStat{2,1}{2,2}{5,2};
    %%Between groups  (e.g.,=f(nb of group))
    %dfbgT1  =GenStat{2,1}{2,2}{4,2};
    %%F-value for the interaction effect
    %FvalueT1=GenStat{2,1}{2,2}{4,5};
    %%p-value of the tested interaction effect
    %pvalueT1=GenStat{2,1}{2,2}{4,6};
    %GenStat2ReportT1={dfbgT1, dfwgT1, FvalueT1, pvalueT1};
    %display(GenStat2ReportT1)
    
    
    %%----------------------------------------------------------------------
    %Stat T1: ANCOVA (effect of FA2 on the slope of data=f(FA3) compared to 1.
    %----------------------------------------------------------------------
    %e.g., FA 1 is coherence
    %e.g., FA 2 is prior's strength
    %e.g., FA 3 is displayed dir
    
    %Loop over each condition FA1-FA2
    %Levels of FA 1
    for k=1:FA.g1.lvlsnb
        %Levels of FA 2
        for j=1:FA.g2.lvlsnb
            %Set x: displayed dirs (levels of FA3)
            x=[];
            x=reshape(repmat(FA.g3.lvlsnm, 2, 1), [],1);
            %Set y: data (e,g., estimated dirs and dirs with slope 1).
            y=[];
            y=reshape([fig1.datamean.deg(j,:,k)', FA.g3.lvlsnm], [],1);
            %Set the groups (g)
            g=[];
            g=reshape(repmat([1,2], FA.g3.lvlsnb,1), [],1);
            %    %Title rows
            %    GenStatAgainst1{j+1,  1}=[FA.g2.nm num2str(FA.g2.lvlsnm(j))];
            %    %title columns
            %    GenStatAgainst1{  1,k+1}=[FA.g1.nm num2str(FA.g1.lvlsnm(k))];
            %run analysis
            Ancova_effect_of_FA2vs1{j+1,k+1}=statAncova(x,y,g);
        end
    end
    
    %Stat to report (write sthg here)
    
    
    %%---------------------------------------------------------------------
    %Stat: ANCOVA - compare the slopes of FA 2's levels for each FA 1.
    %---------------------------------------------------------------------
    %e.g., FA 1 is coherence
    %e.g., FA 2 is prior's strength
    %e.g., FA 3 is displayed dir
    %Loop over the levels of FA1 (e.g., coherences)
    for k=1:FA.g1.lvlsnb
        %Set x: FA3 - e.g., displayed dirs (levels of FA 3)
        x = [];
        x = reshape(repmat(FA.g3.lvlsnm, FA.g2.lvlsnb,1),[],1);
        
        %Set y: data e.g., estimated dirs
        y = [];
        y = reshape(fig1.datamean.deg(:,:,k)', [], 1);
        
        %Set g: groups (levels of FA2, e.g., priors)
        g = [];
        g = reshape(repmat(FA.g2.lvlsnm', FA.g3.lvlsnb,1), [],1);
        
        %title columns
        GenStatBetwLv{1,k} = [FA.g1.nm num2str(FA.g1.lvlsnm(k))];
        
        %run analysis
        GenStatBetwLv{2,k} = statAncova(x,y,g);
    end
    
    %Stat to report(write sthg here)
    
    
    %----------------------------------------------------------------------
    %Store the statistics
    %----------------------------------------------------------------------
    %Title columns
    GenStat=[ {'Ancova_effect_of_FA2'},...
        {'Ancova_effect_of_FA2vs1'},...
        {'Ancova_lvl1_vs_lvl2_FA1'}];
    %Data
    GenStat(2,:) = [ {Ancova_effect_of_FA2},...
        {Ancova_effect_of_FA2vs1},...
        {GenStatBetwLv}];
    fprintf ('%s \n','(slplotMeanAndVarBehavior) The results of the Ancova have been stored')
    
    %Report the statistics of the Ancova
    %Degrees of freedom
    %Within groups  (e.g.,=f(nb of subjects))
    dfwgT1   = GenStat{2,1}{2,2}{5,2};
    %Between g roups  (e.g.,=f(nb of group))
    dfbgT1   = GenStat{2,1}{2,2}{4,2};
    %F-value for the interaction effect
    FvalueT1 = GenStat{2,1}{2,2}{4,5};
    %p-value of the tested interaction effect
    pvalueT1 = GenStat{2,1}{2,2}{4,6};
    GenStatT1 = {dfbgT1, dfwgT1, FvalueT1, pvalueT1};
    display(GenStatT1)
    
    %Report the statistics of the Ancova
    %Degrees of freedom
    %Within groups  (e.g.,=f(nb of subjects))
    dfwgT2  = GenStat{2,2}{2,2}{2,2}{5,2};
    %Between groups  (e.g.,=f(nb of group))
    dfbgT2  = GenStat{2,2}{2,2}{2,2}{4,2};
    %F-value for the interaction effect
    FvalueT2 = GenStat{2,2}{2,2}{2,2}{4,5};
    %p-value of the tested interaction effect
    pvalueT2 = GenStat{2,2}{2,2}{2,2}{4,6};
    GenStatT2 = {dfbgT2, dfwgT2, FvalueT2, pvalueT2};
    display(GenStatT2)
    
    %Report the statistics of the Ancova
    %Degrees of freedom
    %Within groups  (e.g.,=f(nb of subjects))
    dfwgT3  = GenStat{2,3}{2,1}{2,2}{5,2};
    %Between groups  (e.g.,=f(nb of group))
    dfbgT3  = GenStat{2,3}{2,1}{2,2}{4,2};
    %F-value for the interaction effect
    FvalueT3 = GenStat{2,3}{2,1}{2,2}{4,5};
    %p-value of the tested interaction effect
    pvalueT3 = GenStat{2,3}{2,1}{2,2}{4,6};
    GenStatT3 = {dfbgT3, dfwgT3, FvalueT3, pvalueT3};
    display(GenStatT3)
    
end

%Backup figure
autobackup(fig1.hdle,fig1.nm,'.fig');

%[fig]=graphStat(FA, est_dir,fig);
%Graph the quantitative effect of the prior's strength
if FA.g1.lvlsnb ==2
    [fig3] = graphPriorEffect(fig1,FA);
end
