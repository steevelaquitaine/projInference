
%Model's predictions differ the most:
%------------------------------------
%   - When prior strength is ~3xllh strength.
%   - When llh mean is further from the prior mean outside the direction
%   range containing the prior modes and mean.
%
%   - Overlap between prior and llh distributions should be minimal.
%   in that case:
%   - Bayesian posterior tends toward unimodality.
%   - Competition model almost always has 3 readout estimates which get
%clearer when prior and llh modes have the same height (
%prior strength~3xllh strength, look at very low coherence).
%
%note:  llhk is taken from average best fit with competition model with
%llhk for coherence 6%
%     usage: simBayesInferenceVsCompetitionUnimodal([225 215 195 70],[2.2 2.2 2.2 2.2],225,0.86)

function post=simBayesInferenceVsCompetitionUnimodal(llhmean,llhstd,priormode,priorstd)

%%
%likelihood
numCases=length(llhmean);
for i=1:numCases
    llh(:,i)=vmPdfs(1:1:360,llhmean(i),llhstd(i),'norm');
end
%unimodal Gaussian prior
Prior=vmPdfs(1:1:360,priormode,priorstd,'norm');


%Bayesian inference
%------------------
%posterior
for i=1:numCases
    post(:,i)=llh(:,i).*Prior/sum(llh(:,i).*Prior);
end
%draw Bayesian inference
figure('color','w')

for i=1:numCases
    subplot(numCases,2,2*i-1)
    hold all;
    plot(llh(:,i),'k','linesmoothing','on','linewidth',1)
    plot(Prior,'-','linewidth',2,'color',[0.7 0.7 0.7])
    plot(post(:,i),'r','linesmoothing','on','linewidth',2)
    
    %show Bayes MAP readout estimates
    MAPS=find(post(:,i)==max(post(:,i)));
    h(1)=plot(MAPS,post(MAPS,i),'R.','markersize',30);
    
    %Competitions readout estimate
    h(2)=plot(llhmean(i),llh(llhmean(i),i),'.','markersize',30,'color',[.3 .3 1]);
    plot(priormode(1),Prior(priormode(1)),'.','markersize',30,'color',[.3 .3 1])
    %plot(priormode(2),Prior(priormode(2)),'.','markersize',30,'color',[.3 .3 1])
    
    %graphics
    %-------
    xlim([0 360])
    ylabel('Probability')
    xlabel('Motion directions (degrees)')
    
    %title
    if i==1
        title(['Bayesian Inf. & Comp. readout estimates'])
        %legend(h,'Bayes Inf.','Comp')
        legend boxoff
    end 
end





