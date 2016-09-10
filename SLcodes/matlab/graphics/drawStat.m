
      %Author: Steeve Laquitaine
        %Date: 131007
  %last modif: 140418
       %usage:
%             f0=randsample(1:1:360,10000,'true');
%             f1=randsample([0.24 0.12 0.06],10000,'true');
%             f2=randsample([80 40 20 10],10000,'true');
%             data=randsample(1:1:360,10000,'true');
%             [means,stds]=drawStat(data',f0,f1,f2)
%            
%Description:
    %any data
    %f0,f1,f2 are three vectors (factors) associated with the data.

function [means,stds]=drawStat(data,d,coh,pstd,varargin)

%status
fprintf('%12s \n','(drawStat) I am now drawing data statistics...please wait')
tic

%check if only one factor is input, ignore other two factors
if nargin<3
    coh=ones(numel(d),1);
    pstd=ones(numel(d),1);
end

%make sure data is a column vector
if size(data,2)>size(data,1)
    data=data';
end

%make sure d is a row vector
if size(d,1)>size(d,2)
    d=d';
end

%make sure coh is a row vector
if size(coh,1)>size(coh,2)
    coh=coh';
end

%make sure coh is a row vector
if size(pstd,1)>size(pstd,2)
    pstd=pstd';
end


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
c=colormap;
F.f2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};
marksz=70;

%Make mean & std
%mean
figure('color','w')
for j=1:F.f2.n
    h(j)=subplot(1,F.f2.n,j);
    for k=1:F.f1.n
        for i=1:F.f3.n
            means(k,i,j)=nanmean(data(F.inter.pos{k,i,j},:));
        
            %draw
            hold all
            scatter(F.f1.L(k),means(k,i,j)',...
                marksz,...
                'MarkerEdgeColor','w',...
                'MarkerFaceColor',F.f2.color{i},...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
    
            %ideal
            plot(F.f1.L(k),F.f1.L(k),'k.','Markersize',10)
        end
    end
    
    %get axis ylimits
    %case where means are nan
    if ~isnan(max(means(:))) 
        ymax(j) = max(means(:));
    else
        ymax(j) = 1e-10;
    end
    if ~isnan(max(means(:)))
        ymin(j) = min(means(:));
    else
        ymin(j) = 0;
    end
end

%set axis ylimits
set(h,'ylim',[max(ymin) max(ymax)])
% SLremoveDeadSpace(0.05)


%(case we want std too)
if sum(strcmp(varargin,'std'))==1
    %std
    figure('color','w')
    for j=1:F.f2.n
        h(j)=subplot(1,F.f2.n,j);
        for k=1:F.f1.n
            for i=1:F.f3.n
                stds(k,i,j)=nanstd(data(F.inter.pos{k,i,j},:));
                
                %draw
                hold all
                scatter(F.f1.L(k),stds(k,i,j)',...
                    marksz,...
                    'MarkerEdgeColor','w',...
                    'MarkerFaceColor',F.f2.color{i},...
                    'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
                ylabel('Std','fontsize',14)
            end
        end
        
        %get axis ylimits
        ymax(j)=max(stds(:));
        ymin(j)=min(stds(:));
        
        %case where stds are nan
        if ~isnan(max(means(:)))
            ymax(j) = max(stds(:));
        else
            ymax(j) = 1e-10;
        end
        if ~isnan(min(stds(:)))
            ymin(j) = min(stds(:));
        else
            ymin(j) = 0;
        end
    end
    
    %set axis ylimits
    set(h,'ylim',[max(ymin) max(ymax)])
    % SLremoveDeadSpace(0.05)
end


%tell when it ends
fprintf('\n %12s \n','drawing data statistics...done')

toc