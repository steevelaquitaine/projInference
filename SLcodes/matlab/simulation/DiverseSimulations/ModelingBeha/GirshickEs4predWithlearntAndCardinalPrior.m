
 %Author: Steeve Laquitaine
   %date: 140205
%Purpose: draw predictions of the distribution of MAP estimates produced by
%Girshick estimator with cardinal and learnt priors for input 
%cardinal, learnt prior & likelihood strengths.

%Usage
    %k.m=[10 10]; 
    %k.pc=[0 42]; 
    %k.pl=[0 42]; 
    %[mPdfs,l,PRIORcardinal,po,di,m,uniqMAP,whichplot]=GirshickEs4predWithBimodalPrior(k,[0.9 0.65 0.5;
    %0.8 0 0];)
    
%References
%http://mathworld.wolfram.com/vonMisesDistribution.html
%Girshick, A. R., Landy, M. S. & Simoncelli, E. P. Cardinal rules:
%visual orientation perception reflects knowledge of environmental 
%statistics. Nature Publishing Group 14, 926?932 (2011).
%Swindale, N. V. Orientation tuning curves: empirical description 
%and estimation of parameters. Biol Cybern 78, 45?56 (1998).
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/237503

function [mPdfs,l,PRIORcardinal,po,di,m,uniqMAP,whichplot]=GirshickEs4predWithlearntAndCardinalPrior(k,colorM)

%direction(degrees)
%To increase the resolution of MAPs estimates, we can increase the
%motion directions resolution "di" e.g. [1:0.5 360] instead
%of [1 1 360]. It works fine.
di=[1:1:360]';

%measurements(degrees)
m=[1:1:360];
n=numel(m);

%measurement densities with mean di (col) over measurement m(row)
mPdfs=vmPdfs(m,di,k.m,'norm');

%likelihood of measurement mi(col) given motion directions di (row)
ul=[1:1:360];
l=vmPdfs(di,ul,k.m,'norm');

%cardinal PRIOR (over motion directions (row) is same for each mi(col))
PRIORcardinal=vmPdfs(di,[90 180 270 360],k.pc,'norm');
PRIORcardinal=0.25.*sum(PRIORcardinal,2);
PRIORcardinal=PRIORcardinal(:,ones(numel(m),1));

%learnt PRIOR (over motion directions (row) is same for each mi(col))
up=225;
PRIORlearnt=vmPdfs(di,up,k.pl,'norm');
PRIORlearnt=PRIORlearnt(:,ones(numel(m),1));

%get posteriors the common cause density of the 3 representations.
%avoid computing errors at <10^10 floating points;
po=PRIORcardinal.*PRIORlearnt.*l;
Zpo=sum(po,1);
Zpo=Zpo(ones(numel(di),1),:);
po=po./Zpo;
po=fix(po.*10^10)/10^10;

%find MAPs estimates the most likely common causes (values) of each ul/mi 
%(rows). Same mi may produce more than one map. e,g., at the obliques 
%(col). And different mi may produce same map.
MAP=nan(numel(m),numel(di));
for i=1:numel(m)
    numMAPs=sum(po(:,i)==max(po(:,i)));
    MAP(i,1:numMAPs)=di(po(:,i)==max(po(:,i)));
end

%get the max likelihood of observing each data "upo" that is the probability 
%of ul/mi given the motion direction "di" for each "ul" and associaed "upo". 
%ul/mi
maxnumMAP=max(sum(~isnan(MAP),2));
ul=m';
ul=ul(:,ones(maxnumMAP,1));
ul=ul(:);

%associated MAP estimate
MAP=MAP(:,1:maxnumMAP);
MAP=MAP(:);

%associated likelihood for each motion direction (column:1:1:360 degrees)
likelihoodMAPgivenMi=repmat(mPdfs,maxnumMAP,1);

%everything sorted by MAP
TheMatrix=[ul MAP likelihoodMAPgivenMi];
TheMatrix=sortrows(TheMatrix,2);
ul=TheMatrix(:,1);
MAP=TheMatrix(:,2);
likelihoodMAPgivenMi=TheMatrix(:,3:end);

%get the likelihood of each MAP (if a MAP had been produced by different 
%mi, then its likelihood is the average probability of this mi for each
%motion direction, it is consistent with the laws of probability). It is
%the (1/nummi)*P(MAP/mi1) + (1/nummi)*P(MAP/mi2) +...
%It works very nice, the estimates density obtained are very smoothed.
uniqMAP=unique(MAP(~isnan(MAP)));
likelihoodUniqMAPgivenMi=nan(numel(uniqMAP),size(likelihoodMAPgivenMi,2));
for i=1:numel(uniqMAP)
    posUniqMAP=find(MAP==uniqMAP(i));
    likelihoodUniqMAPgivenMi(i,:)=mean(likelihoodMAPgivenMi(posUniqMAP,:),1); 
end

%Set the likelihood of MAPs not produced at 0 because the model cannot
%produce those estimates even at a reasonable resolutions of motion 
%direction
MAPnotProduced=setdiff(1:1:360,uniqMAP)';
likelihoodMAPnotProduced=zeros(numel(MAPnotProduced),size(mPdfs,2));
      
%Add not produced MAP and their 0 likelihood & sort everything again by MAP
%We do not plot the 0 likelihood of the non produced maps for visualization
%purpose.
likelihoodMAPs=[likelihoodUniqMAPgivenMi;likelihoodMAPnotProduced];
allMAP=[uniqMAP;MAPnotProduced];
TheMatrixII=[allMAP likelihoodMAPs];
TheMatrixII=sortrows(TheMatrixII,1);
allMAP=TheMatrixII(:,1);
likelihoodMAPs=TheMatrixII(:,2:end);


%Plots
whichplot=drawEstimatePdfs(di,m,uniqMAP,mPdfs,likelihoodUniqMAPgivenMi,colorM);

function whichplot=drawEstimatePdfs(di,m,uniqMAP,mPdfs,likelihoodUniqMAPgivenMi,colorM)

%graphics
markSz=8;

%Draw measurement densities and their associated MAP estimates
%densities for many example motion directions.
motion=[45 55 65 75 90];
for i=1:numel(motion)
    motionpos=find(di==motion(i));
    
    %draw measurement density
    subplot(2,numel(motion),i)
    hold all
    title(strcat(num2str(motion(i)),' degrees motion'),'fontsize',14,...
        'fontweight','Bold')
    plot(m,mPdfs(:,motionpos),'.','color',[.7 .7 .7],...
        'linesmoothing','on',...
        'Markersize',markSz)
    maxP=max(mPdfs(:,motionpos));
    %motion direction
    plot([motion(i) motion(i)],[0 maxP],'--','color',[.4 .4 .4],'linewidth',2)
    xlabel('Measurements (degrees)','fontsize',14)
    ylabel('Probability','fontsize',14)
    xlim([1 360])
    set(gca,'xtick',[1 90 180 270 360],'xticklabel',[1 90 180 270 360],'fontsize',14)
    box off
    
    %draw MAP estimate density
    subplot(2,numel(motion),i+numel(motion))
    hold all
    whichplot=plot(uniqMAP,likelihoodUniqMAPgivenMi(:,motionpos),'.','color',colorM,...
        'linesmoothing','on',...
        'Markersize',markSz);
    
    %motion direction
    maxP=max(likelihoodUniqMAPgivenMi(:,motionpos));
    plot([motion(i) motion(i)],[0 maxP],'--','color',[.4 .4 .4],'linewidth',2)
    xlabel('MAP estimates (degrees)','fontsize',14)
    ylabel('Probability','fontsize',14)
    xlim([1 360])
    set(gca,'xtick',[1 90 180 270 360],'xticklabel',[1 90 180 270 360],'fontsize',14)
    box off
    
end
