% slsimulateBayesianModelEstimateSpace.m
%
%     author: steeve laquitaine
%       date: 140522 (last update 141121)
%    purpose: 
%      usage:
%
%
%       slsimulateBayesianModelEstimateSpace(5:20:360)


function slsimulateBayesianModelEstimateSpace(evidence)

%evidence
subplot(311)
title('Evidence')
colors = linspecer(length(evidence));
for i = 1 : length(evidence)
    hold all
    plot(evidence(i),zeros(1,length(evidence)),'.','markersize',30,'color',colors(i,:))
end
xlim([0 360])

%predicted circular estimate space
subplot(312)
title({'Posterior and estimates in circular space: Bias is weaker away',...
    'from prior mean and variability larger as evidence sampled from',... 
    'estimate distributions more often fall either clockwise or', ...
    'counterclockwise to the prior mean which leads their associated estimates',...
    'to more often lie either clockwise or counterclockwise to the prior mean',... 
    'reducing bias and increasing estimate dispersion. This effect is due to the circularity of the space.'})
slsimulateBayesianModelCircularEstimateSpace(5:20:360,5:10:360,5,225,4.77,0,'vonMisesPrior','withoutCardinal');
ylim([0 0.22])

%predicted linear estimate space
subplot(313)
title('Posterior and estimates in linear space: The effect disappears on linear space',...
    'predicting constant variability with increasing distance to prior mean')
slsimulateBayesianModelLinearEstimateSpace(5:20:360,5:10:360,30,225,30,0,'vonMisesPrior','withoutCardinal');
ylim([0 0.22])


%add plots of variability