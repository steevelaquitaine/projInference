

%usage:
%       
%     tau = [0.29 0.08 4.00 18.06 20.00 20.35 1.68 20.00 21.20 55.11 16.69 1.57];
%     slPlotBayesianLearningRates(tau,200)


function slPlotBayesianLearningRates(tau,nTrials)

%Trials
t = 1 : nTrials;

%colors
colr =  linspecer(length(tau));
figure('color','w')
box off
for i = 1 : length(tau)
    hold on; 
    plot(1-exp(-t/tau(i)),'color',colr(i,:))
    ylim([0 1])
end
ylabel('Prior strength learning (% of final strength)')
xlabel('Trials (%)')