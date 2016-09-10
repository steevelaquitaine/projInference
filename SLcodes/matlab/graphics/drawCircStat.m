
%SLdrawCircStat.m
%
%     Author: Steeve Laquitaine
%       Date: 131007 last modif: 140601
%
%      usage:
%        f0=randsample(1:1:360,1000,'true');
%        f1=randsample([0.24 0.12 0.06],1000,'true');
%        f2=randsample([80 40 20 10],1000,'true');
%        data=f0;
% 
%        [means,stds]=drawCircStat(data,f0,f1,f2)
% 
%
%
% Description:
%   data are in degree.
%   f0,f1,f2 are three vectors of conditions associated with the data.
% 
% noteTODO: 
%   add an "errorbar" option where data can be plot like with
%   errorbar in the plot for the means instead of with just dots (figure 1).

function [means,stds] = drawCircStat(data,d,coh,pstd,plotType,varargin)


%If only one factor is input, ignore other two factors
if isempty(coh)
    coh=ones(numel(d),1);
end
if isempty(pstd)
    pstd=ones(numel(d),1);
end

%data(cartesians)
data=polar2cartesian(data,1);

%factors 1,2,3
F.f1.i=d;
F.f1.nm='d';
F.f1.L=unique(F.f1.i);
F.f1.L=sort(F.f1.L,'ascend');
F.f1.n=numel(F.f1.L);

F.f2.i=coh;
F.f2.nm='coh';
F.f2.L=unique(F.f2.i);
F.f2.L=sort(F.f2.L,'descend');
F.f2.n=numel(F.f2.L);

F.f3.i=pstd;
F.f3.nm='pstd';
F.f3.L=unique(F.f3.i);
F.f3.L=sort(F.f3.L,'descend');
F.f3.n=numel(F.f3.L);

%positions main
for i=1:F.f1.n
    F.f1.pos(i)={find(F.f1.i==F.f1.L(i))};
end
for i=1:F.f2.n
    F.f2.pos(i)={find(F.f2.i==F.f2.L(i))};
end
for i=1:F.f3.n
    F.f3.pos(i)={find(F.f3.i==F.f3.L(i))};
end

%positions inter
for k=1:F.f1.n
    for j=1:F.f2.n
        for i=1:F.f3.n
            F.inter.pos(k,i,j)=...
                {intersect( ...
                intersect(F.f1.pos{k},F.f2.pos{j}),...
                F.f3.pos{i})};
        end
    end
end

%Graphics
F.f2.color={[0.5 0 0],...
    [1 0 0],...
    [1 0.5 0],...
    [0.65 0.65 0]};
marksz=100;

%---------
%Make mean
%---------
if ~exist('plotType','var') 
    figure('color','w')
end
means=[];
stds=[];
for j=1:F.f2.n
    if ~exist('plotType','var') 
        subplot(1,F.f2.n,j)
        axis square
    end
    for k=1:F.f1.n
        for i=1:F.f3.n
            
            %get mean
            stat{k,i,j} = SLcircMeanStd(data(F.inter.pos{k,i,j},:));
            means(k,i,j) = stat{k,i,j}.deg.mean;
            
            %plot without errorbars
            if ~exist('plotType','var')
                hold all
                scatter(F.f1.L(k),means(k,i,j)',...
                    marksz,...
                    'MarkerEdgeColor','w',...
                    'MarkerFaceColor',F.f2.color{i},...
                    'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
            end
            
            if ~exist('plotType','var')
                
                %ideal
                plot(F.f1.L,F.f1.L,'k--','Markersize',20,'linewidth',1)
                ylabel('Mean (degrees)','fontsize',14)
                axis square
                
                %graphics
                %--------
                set(gca,...
                    'xtick', 0:45:360,'xticklabel',0:45:360,...
                    'ytick', 0:45:360,'yticklabel',0:45:360,...
                    'fontsize',13);
                xlim([0 360]);
                ylim([0 360]);                
            end
        end
    end
    xlim([0 360])
    ylim([0 360])
    
    %make clean
    drawPublishAxis
end

if ~isempty(varargin) && strcmp(varargin{:},'clean=on')
    SLremoveDeadSpace
end



%---
%std
%---
h=[];
maxplot=[];
if ~exist('plotType','var') 
    figure('color','w')
    maxplot=[];
end
for j=1:F.f2.n
    if ~exist('plotType','var') 
        subplot(1,F.f2.n,j)
        axis square
    end
    for k=1:F.f1.n
        for i=1:F.f3.n
            
            %get std
            stds(k,i,j)=stat{k,i,j}.deg.std;
            
            %plot without errorbars
            if ~exist('plotType','var')
                hold all
                scatter(F.f1.L(k),stds(k,i,j)',...
                    marksz,...
                    'MarkerEdgeColor','w',...
                    'MarkerFaceColor',F.f2.color{i},...
                    'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
                ylabel('Std (degrees)','fontsize',14)

                %graphics
                %--------
                %x-axis
                set(gca,'xtick',0:45:360,'xticklabel',0:45:360,...
                    'fontsize',13)
                set(gca,'ytick',0:10:100,'yticklabel',0:10:100,'fontsize',13);
                xlim([0 360]);
                
                %get max plot
                h=[h gca];
                maxplot=[maxplot max(stds(k,i,j))];
            end
        end
    end
    
    %plot without errorbars
    if ~exist('plotType','var')
        xlim([0 360])
        ylim([0 360])
    end
    
end
for i=1:length(h)
    ylim(h(i),[0 max(maxplot)])
end
if ~isempty(varargin) && strcmp(varargin{:},'clean=on')
    SLremoveDeadSpace
end

%---------------------
%Or plot with errorbars
%----------------------
if exist('plotType','var') && strcmp(plotType,'errorbar')
    figure('color','w')
    h=[];
    maxplot=[];
    box off
    for j=1:F.f2.n
        subplot(1,F.f2.n,j)
        axis square
        for k=1:F.f1.n
            for i=1:F.f3.n
                myerrorbar(F.f1.L(k),means(k,i,j)',...
                    'yError',stds(k,i,j),...
                    'Symbol=o',...
                    ['MarkerFaceColor=',num2str(F.f2.color{i})])
            end
            
            %ideal
            plot(F.f1.L,F.f1.L,'k--','Markersize',20,'linewidth',1)
            ylabel('Mean (degrees)','fontsize',14)
            axis square
            set(gca,'xtick',0:30:360,'xticklabel',0:30:360)
            set(gca,'ytick',0:30:360,'yticklabel',0:30:360)
            set(gca,'fontsize',14)
            xlim([0 360])
            ylim([0 360])      
            
            %get max plot
            h=[h gca];
            %maxplot=[maxplot max(stds(k,i,j))];
        end
    end
end
%for i = 1 : length(h)
    %ylim(h(i),[0 max(maxplot)])
%end
if ~isempty(varargin) && strcmp(varargin{:},'clean=on')
    SLremoveDeadSpace
end
