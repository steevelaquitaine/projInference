%SLmakeDicreteCountVMdensity.m
%
%  author: steeve laquitaine
%    date: 2015/05/06
% purpose: create a discrete von Mises count histogram (y axis are integers)
%
%   usage:
%
%        [series,count] = SLmakeDicreteCountVMdensity(45:60:355,225,0.74848,143)
 
function [series,count] = SLmakeDicreteCountVMdensity(mySamples,myMean,K,nTrials)
%We make the probability densities we need based on von mises or gaussian
%laws. Then we calculate the number of occurrence of each direction
%(pdf*numTrials). Occurences are not real numbers because pdf are not but
%we need real numbers. We could randomly sample the density but
%for small sample size as it is the case in our experiment, this method
%creates small left/right asymetries relative to the mean. Those asymetries
%could create confounding biases in estimation that are not accounted for
%by the true continuous prior density.
%We avoid this asymetry by rounding the number of occurrence to real values
%and then by replicating direction trials according to the resulting number
%of occurrence. The incovenient is that the parameters of the prior std
%change but only slightly.
series = [];
count  = [];
nSamples = length(mySamples);

%PriorShape = gauss_distribution(mySamples.degree,mean,K);
PriorShape = vmPdfs(mySamples,myMean,K,'norm');

%create trials' count. We repeat samples according to the density.
contDistInCount = PriorShape*nTrials ;

for i = 1 : numel(mySamples)
    
    contDistInCountrepet = repmat(mySamples(i),round(contDistInCount(i)),1);
    series = [series; contDistInCountrepet];
    
    clear contDistInCountrepet
    
end
series = series';


%occurrence of displayed feature
for i = 1 : nSamples
    count(i) = numel(find(series==mySamples(i)));
end
%
% %or generate a uniform prior
% if K==inf
%     %check if "nTrials" is a multiple of
%     %"nSamples".
%     if rem(nTrials,nSamples)==0
%         fprintf(['--- A UNIFORM distribution is being drawn ----'])
%         count=repmat(nTrials/nSamples,...
%             1,nSamples);
%         series=repmat(mySamples,count(1),1);
%         series=series(:);
%         series=series';
%     else
%         fprintf(['\n (SLinitRunUniPriorLocTask) "nTrials" must be a multiple of "nSamples" \n'])
%         return
%     end
% end