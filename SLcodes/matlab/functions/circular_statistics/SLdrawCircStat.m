
%SLdrawCircStat.m
%
%
% Author: Steeve Laquitaine
%   Date: 131007
%purpose: draw a circular data mean and std against circular
%         independent variables "v1" and up to two more variables "v2","v3"
%         (e.g., dataRaw = mean motion direction estimates
%         and v1 = true motion direction, v2 = motion coherence, v3 = 
%         direction statistical distribution strength)
%
% usage:
%
%     v1 = randsample(1:1:360,1000,'true');
%     v2 = randsample([0.24 0.12 0.06],100,'true');
%     v3 = randsample([80 40 20 10],1000,'true');
%     data = v1;
%     [means,stds] = SLdrawCircStat(data,v1,v2,v3)
%
%     %or in two subplots  
%     [means,stds] = SLdrawCircStat(data,v1,v2,v3,'mysubplots','subplot(1,2,1)')
%     [means,stds] = SLdrawCircStat(data,v1,v2,v3,'mysubplots','subplot(1,2,2)')
%
% Outputs : 
%
%       means : means by condition
%       stds : std by condition
%       linfit : 
%           linfit.slope_fitMean : slope from fitting tthe raw data
%           linfit.interc_fitMean : intercept
%           linfit.slope_boot : boostrapped slopes
%           linfit.interc_boot : bootstrapped intercepts
%           linfit.coh : coherences for each slope and intercept
%           linfit.prior : priors 
%       stat_boot : 
%        feat_lin : v1 linearized for fitting and plot
%          est_lin: means linearized for fitting and plot
%
% Description:
%     - data are in degree.
%     - data mean "y" are linearized for the scatter plot  as linear distance
%     - from "x" by "y = x + signed angle(y,x)"
%     - linear fits are performed on linearized trial-data in the same way
%
%       I linearize because Linear fit of raw data is not possible because data
%       are circular and certain data will look like outliers in the linear
%       space (e.g., 1 deg estimate for 360 deg motion direction)
%       but are not in the circular space.
%
%varargin:
%
%       e.g., 'mysubplots', 'subplot(1,1,1)': set your subplots externally
%                                             otherwise default
%       e.g., 'mycolor',[1 0 0]: set your color
%       e.g., 'bootstrapSlopeAndIntercept': bootstrap linear fit slope and
%                                           intercept

%draw circular statistics
function [means,stds,linfit,stat_boot,feat_lin,est_lin] = SLdrawCircStat(dataRaw,v1,v2,v3,varargin)

%initialize
stat_boot = [];

%data(cartesians)
data = SLpolar2cartesian(dataRaw,1);

%factors 1,2,3
F.f1.i = v1;
F.f1.nm = 'v1';
F.f1.L = unique(F.f1.i);
F.f1.L = sort(F.f1.L,'ascend');
F.f1.n = numel(F.f1.L);

%case data are plotted against two variables
if ~isempty(v2)
    F.f2.i = v2;
    F.f2.nm = 'v2';
    F.f2.L = unique(F.f2.i);
    F.f2.L = sort(F.f2.L,'descend');
    F.f2.n = numel(F.f2.L);
else
    F.f2.i = ones(size(v1));
    F.f2.nm = 'None';
    F.f2.L = 1;
    F.f2.n = 1;
    F.f3 = F.f2;
end

%case data are plotted against three variables
if ~isempty(v3)
    F.f3.i = v3;
    F.f3.nm='v3';
    F.f3.L=unique(F.f3.i);
    F.f3.L=sort(F.f3.L,'descend');
    F.f3.n=numel(F.f3.L);
else
    F.f3.i = ones(size(v1));
    F.f3.nm = 'None';
    F.f3.L = 1;
    F.f3.n = 1;
end

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
c = colormap;

%set default subplot or input subplot
F.f2.color = {[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};
if any(strcmp(varargin,'mycolor'))
    F.f2.color = varargin(find(strcmp(varargin,'mycolor'))+1);
    szr = size(F.f2.color{:});
    %case one color
    if szr(1) == 1
        F.f2.color = repmat(F.f2.color,4,1);
    end
end

%====
%mean
%====
figure(1)
set(gcf,'color','w')
for j = 1 : F.f2.n
    %set default subplot or input subplot
    if ~any(strcmp(varargin,'mysubplots'))
        subplot(1,F.f2.n,j)
    else
        eval(varargin{find(strcmp(varargin,'mysubplots'))+1})
    end
    %visualize data
    for i = 1 : F.f3.n
        for k = 1 : F.f1.n
            stat{k,i,j} = SLcircMeanStd(data(F.inter.pos{k,i,j},:),'cartesian');
            means(k,i,j) = stat{k,i,j}.deg.mean;
            %data
            %linearize data means for scatter plot (otherwise 1 deg looks
            %further from 360 deg). The max angular distance between the y
            %and x are within -180 and 180 degrees.
            [feat_lin(k),est_lin(k,i,j)] = slLinearizeCircDataForScatterAndLinFit(F.f1.L(k),means(k,i,j));
            
            %plot errorbar
            hold on;
            myerrorbar(feat_lin(k),est_lin(k,i,j)','yError',stat{k,i,j}.deg.std,...
                'MarkerEdgeColor=w',...
                ['MarkerFaceColor=' num2str(F.f2.color{i})]);
            
        end
        %ideal prediction
        plot(1:1:360,1:1:360,'k:','Markersize',3)
        
        %linearize trial data for linear fit of trial data
        est = dataRaw(F.f3.i==F.f3.L(i) & F.f2.i==F.f2.L(j));
        feat = v1(F.f3.i==F.f3.L(i) & F.f2.i==F.f2.L(j));
        [feat,est] = slLinearizeCircDataForScatterAndLinFit(feat,est);
        
        %=====================================
        %get slope, intercept fitted mean data
        %=====================================
        [a,~] = sllinefit(feat,est,1:1:360,F.f2.color{j});
        linfit.slope_fitMean(j,i) = a(1);
        linfit.interc_fitMean(j,i) = a(2);
        
        %=============================
        %bootstrap slope and intercept
        %=============================
        if any(strcmp(varargin,'bootstrapSlopeAndIntercept'))                        
            x=feat; y=est;displayOp='''dispOff''';                       
            myfun = ['stat=sllinefit(x,y,1:1:360,[' num2str(F.f2.color{j}) '],' displayOp ');'];
            stat_boot = slBootstrap(myfun,x,y,length(x),1000);
            stat_boot = cell2mat(stat_boot');
            linfit.slope_boot = stat_boot(:,1);
            linfit.interc_boot = stat_boot(:,2);
        end
    end
    xlim([0 360])
    ylim([0-180 360+180])
end

%===
%std
%===
figure(2)
set(gcf,'color','w')
for j = 1 : F.f2.n
    %set default subplot or input subplot
    if ~any(strcmp(varargin,'mysubplots'))
        subplot(1,F.f2.n,j)
    else
        %input subplot
        eval(varargin{find(strcmp(varargin,'mysubplots'))+1})
    end
    %visualize data
    for k = 1 : F.f1.n
        for i = 1 : F.f3.n
            stds(k,i,j)=stat{k,i,j}.deg.std;
            %data
            hold all
            scatter(F.f1.L(k),stds(k,i,j)',...
                'MarkerEdgeColor','w',...
                'MarkerFaceColor',F.f2.color{i},...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
        end
    end
    xlim([0 360])
    ylim([0 100])
end