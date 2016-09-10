
  %Date: 131007
%Author: Steeve Laquitaine
%usage:
%     f0=randsample(1:1:360,10000,'true');
%     f1=randsample([0.24 0.12 0.06],10000,'true');
%     f2=randsample([80 40 20 10],10000,'true');
%     data=f0;
%     [means,stds]=drawCircStat(data,f0,f1,f2)

%Description
    %data is in degree.

%draw circular statistics
function [means,stds]=drawCircStat(data,d,coh,pstd)
%Inputs: 
    %data
    %3 vectors for 3 factors. Rows are level values.

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
        
%Make mean & std
c=colormap;
F.f2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};


%mean
figure('color','w')
for j=1:F.f2.n
    subplot(1,3,j)
    for k=1:F.f1.n
        for i=1:F.f3.n
            stat{k,i,j}=circMeanStd(data(F.inter.pos{k,i,j},:));
            means(k,i,j)=stat{k,i,j}.deg.mean;
        
            %draw
            hold all
            scatter(F.f1.L(k),means(k,i,j)',...
                'MarkerEdgeColor','w',...
                'MarkerFaceColor',F.f2.color{i},...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
            %ideal
            plot(F.f1.L(k),F.f1.L(k),'k.','Markersize',3)
        end
    end
    xlim([0 360])
    ylim([0 360])
end


%std
figure('color','w')
for j=1:F.f2.n
    subplot(1,3,j)
    for k=1:F.f1.n
        for i=1:F.f3.n
            stds(k,i,j)=stat{k,i,j}.deg.std;       
        
            %draw
            hold all
            scatter(F.f1.L(k),stds(k,i,j)',...
                'MarkerEdgeColor','w',...
                'MarkerFaceColor',F.f2.color{i},...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))

        end
    end
    xlim([0 360])
    ylim([0 360])
end



toc