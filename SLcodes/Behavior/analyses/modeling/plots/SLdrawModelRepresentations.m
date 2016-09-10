
% SLdrawModelRepresentations.m
%
%     author: steeve laquitaine
%       date: 140601
%
%    purpose: draw Bayesian model's representations
%
%
% The model have 9 fit parameters
% - 3 von Mises llh concentration parameters k (0<kl<inf, stimulus noise)
% - 4 von Mises (or mixture of von Mises) priors concentration parameters 
%   (0<kp<inf)
% - 1 parameter for the fraction of trial with random estimation (0<prand<1)
% - 1 parameter for motor noise (0<km<inf)
%
%The model can use:
% - von Mises priors
% - bimodal priors (mixture of von Mises priors)
%
%Additionally, the model can accounts for: 
% - cardinal biases (additional Bayesian inference with cardinal priors)
% - random estimation
% - motor noise
%
%
%
%references: 
%     -Hurliman et al, 2002,VR
%     -Stocker&Simoncelli,2006,NN
%     -Girshick&Simoncelli,2011,NN
%     -Chalk&Series,2012,JoV
%
%   usage: SLdrawModelRepresentations(coh,pstd,priorModes,d,fitP)

function output=SLdrawModelRepresentations(coh,pstd,priorModes,d,fitP,PestimateGivenModel,TheModel,varargin)

%get fit parameters. They should be in this order.
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kcardinal','Prand','km'
kcardinal=fitP(8);
% Prandom=simP(9);
% km=simP(10);

%(case we don't want to fit the cardinal prior)
%----------------------------------------------
if strcmp(TheModel,'withoutCardinal')
    kcardinal=0;
end    

%set strength of trial-representations
%likelihood strength
klthisT(coh==0.24)=fitP(1);
klthisT(coh==.12)=fitP(2);
klthisT(coh==.06)=fitP(3);

%prior strength
KlearntthisT(pstd==80)=fitP(4);
KlearntthisT(pstd==40)=fitP(5);
KlearntthisT(pstd==20)=fitP(6);
KlearntthisT(pstd==10)=fitP(7);

%---------------------
%Draw representations
%---------------------
%%% Weakest prior -  llh gets weaker
%-----------------------------------
figure('color','w')
clf
for idx=1:3
    
    %axis
    output.h(idx)=subplot(2,4,idx);
    hold all
    box off
    
    %most distant direction for this prior
    if sum(strcmp(varargin{1},'vonMisesPrior'))==1
        dtoSee = 135;
    end
    if  sum(strcmp(varargin{1},'bimodalPrior'))==1   
        dtoSee=min(d(KlearntthisT==fitP(4)));
    end
    
    %each condition
    j=find(d==dtoSee & KlearntthisT'==fitP(4) & klthisT'==fitP(idx));
    j=j(1);
    i=j;
    
    %learnt prior
    prior=vmPdfs(1:1:360,225,KlearntthisT(i),'norm');
    %plot(prior,'color',[0 0.5 1],'linestyle','-','linewidth',1);
    area(prior,'facecolor',[0 0.5 1],'linestyle','-','linewidth',1);

    %cardinal prior (over motion directions (row) is same for each mi(col))
    PRIORcardinal=vmPdfs(1:1:360,[90 180 270 360],kcardinal,'norm');
    PRIORcardinal=0.25.*sum(PRIORcardinal,2);
    %plot(PRIORcardinal,'color',[.7 .7 .7],'linestyle','-','linewidth',1);
    area(PRIORcardinal,'facecolor',[.7 .7 .7],'linestyle','-','linewidth',1);

    %llh
    llh=vmPdfs(1:1:360,d(i),klthisT(i),'norm');
    %plot(llh,'color',[.5 .5 .5],'linesmoothing','on');
    area(llh,'facecolor',[.5 .5 .5]);

    %posterior
    post=prior.*llh.*PRIORcardinal/sum(prior.*llh.*PRIORcardinal);
    %plot(post,'g','linesmoothing','on');
    area(post,'facecolor','g');

    %estimate density predictions
    %plot(PestimateGivenModel(:,i),'r','linesmoothing','on','linewidth',2);
    area(PestimateGivenModel(:,i),'facecolor','r','linestyle','-','linewidth',1);
    
    %mean and std estimate density
    estimate=(1:1:360)';
    data = SLcircWeightedMeanStd(estimate,PestimateGivenModel(:,i));
    meanEstDensity(i)= data.deg.mean; 
    stdEstDensity(i) = data.deg.std;
        
    %scale plots
    maxP(idx)=max([prior;llh;post;PestimateGivenModel(:,i)]);
    
    %labels
    if idx==1
        ylabel('Probability')
    end
    
    %title and infor
    title({['kl:',num2str(fitP(idx)),'- kp:',num2str(fitP(4))],...
        ['p(es) mean:',num2str(fix(meanEstDensity(i))),',',...
         ' std:',num2str(fix(stdEstDensity(i)))]})

end
for idx=1:3
    ylim(output.h(idx),[0 max(maxP)])
end

%%% weakest llh - strengthening prior
%-----------------------------------
for idx=1:4
    h2(idx)=subplot(2,4,idx+4);
    idxpr=3+idx;
    title(['kl:',num2str(fitP(3)),'- kp:',num2str(fitP(idxpr))])
    hold all
    box off
    
    %most distant direction for this prior
    dtoSee=min(d(KlearntthisT==fitP(idxpr) ));
    
    %get condition
    j=find(d==dtoSee & KlearntthisT'==fitP(idxpr) & klthisT'==fitP(3));
    j=j(1);
    i=j;
    
    %learnt prior
    prior=vmPdfs(1:1:360,225,KlearntthisT(i),'norm');
    plot(prior,'color',[0 0.5 1],'linestyle',':','linewidth',2);
    
    %cardinal prior (over motion directions (row) is same for each mi(col))
    PRIORcardinal=vmPdfs(1:1:360,[90 180 270 360],kcardinal,'norm');
    PRIORcardinal=0.25.*sum(PRIORcardinal,2);
    plot(PRIORcardinal,'color',[.7 .7 .7],'linestyle',':','linewidth',2);
    
    %llh
    llh=vmPdfs(1:1:360,d(i),klthisT(i),'norm');
    plot(llh,'color',[.5 .5 .5],'linesmoothing','on');
    
    %posterior
    post=prior.*llh.*PRIORcardinal/sum(prior.*llh.*PRIORcardinal);
    plot(post,'r','linesmoothing','on');
    
    %estimate density predictions
    plot(PestimateGivenModel(:,i),'r','linesmoothing','on','linewidth',2);
    
    %scale plots
    maxP(idx)=max([prior;llh;post;PestimateGivenModel(:,i)]);
     
    %labels
    if idx==1
        ylabel('Probability')
    end
    if idx==3
        xlabel('motion direction')
    end
end
for idx=1:4
    ylim(h2(idx),[0 max(maxP)])
end
legend('learnt prior','card prior','llh','estimate density') 

