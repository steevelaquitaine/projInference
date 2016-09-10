%modelComparison.m
%
%       $Id: modelComparison.m 750 2012-09-08 07:12:46Z steeve $
%     usage: modelComparison
%        by: steeve laquitaine
%      date: 140521 
%   purpose: compare quantitatively Bayesian inference and Competition 
%prediction for a bimodal prior.
%
%
%     usage:
% meanloglModel1S1=-4.62608263929087;
% meanloglModel2S1=-4.55327094837954;
% meanloglModel1S2=-4.8;
% meanloglModel2S2=-4.2;
% fitPModel1=9;
% fitPModel2=9;
% modelComparison({'sub01','sub02'},[meanloglModel1S1 meanloglModel2S1;...
% meanloglModel1S2 meanloglModel2S2],...
% [fitPModel1 fitPModel2],...
% {'Bayes','Competition'});

function SLmodelComparison(subjects,meanlogl,fitP,modelnm)

%AIC for both models.
%--------------------
%"AIC=2k - 2ln(L)" 
%note that we use expected loglikelihood over trials instead of likelihood
%that we couldn't evaluate because the factorial product of so many 
%probabilities gives 0 as a likelihood.
%AICc is better to use but because it requires likelihood, we cannot use 
%calculate it.
numModels=length(fitP);
numSubjects=length(subjects);
for i=1:length(subjects)
    for j=1:numModels
        AIC(i,j)=2*fitP(1)-2*meanlogl(i,j);
    end
end

%draw AICs for the different models for each subject
figure('color','w')
colors=[0.2 0.2 0.2;
    .7 .7 .7];

%position data
x1(1,:)=1:numModels+1:numModels*numSubjects;
x1(2,:)=x1(1,:)+1;

for i=1:numSubjects
    for j=1:numModels
    hold all
    
    %draw AIC
    b(j)=bar(x1(j,i),AIC(i,j),...
        'edgecolor','none',...
        'facecolor',colors(j,:),...
        'displayname',modelnm{j});
    
    %annotate
    text(x1(j,i),1.01*AIC(i,j),num2str(fix(AIC(i,j)*100)/100),...
        'HorizontalAlignment','center','fontsize',14)
    end
end

%legend
ylim([.9*max(AIC(:)) 1.01*max(AIC(:))])
ylabel('AIC','fontsize',14)
set(gca,'fontsize',14)
xlabel('Subjects','fontsize',14,'fontweight','Bold')

%indicate which subjects
set(gca,'xtick',x1(1,:),'xticklabels',subjects)

%legend
legend(b,'Location','BestOutside')
legend('boxoff')

