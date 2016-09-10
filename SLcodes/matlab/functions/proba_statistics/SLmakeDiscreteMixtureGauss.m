
%SLmakeDiscreteMixtureGauss.m
%
% author: Steeve Laquitaine
%   date: 150109
%purpose: Create a time series of directional data drawn from a discrete
%mixture of Gaussian
%
%usage:
%
%       task = SLmakeDiscreteMixtureGauss([8.4 8.4],5:10:355,[105 345],107)
  
function o = SLmakeDiscreteMixtureGauss(strength,sampleDeg,modes,trialnum)

%initialize
series = [];
count  = [];
sampsiz = numel(sampleDeg);

%Bimodal prior
if strength~=0
    
    %Mixture of two von Mises
    MoG = SLmixtureGaussian(sampleDeg,modes,strength); 
    
    %count of each direction
    f2 = MoG * trialnum ;
    
    %Repeat and store samples
    for i = 1 : numel(sampleDeg)
        f2repet = [];
        f2repet = repmat(sampleDeg(i),round(f2(i)),1);
        series = [series; f2repet];
    end
    series = series';
end

%count
for i = 1 : sampsiz    
    count(i) = numel(find(series==sampleDeg(i)));
end
    
%Uniform distribution
if strength==inf
    
    %"trialnum" must be a multiple of "sampsiz".
    if rem(trialnum,sampsiz)==0
        sprintf('(SLmakeDiscreteMixtureGauss)--- A UNIFORM distribution is being drawn ----')
        count = repmat(trialnum/sampsiz,...
            1,sampsiz);
        series = repmat(sampleDeg,count(1),1);
        series = series(:);
        series = series';
    else
        sprintf('(SLmakeDiscreteMixtureGauss) "trialnum" must be a multiple of "sampsiz"')
        keyboard
    end
end
CirCmean = SLcircMeanStd(SLmakeColumn(series));

%output
o.parameter.dir.sample.degree = sampleDeg;
o.parameter.dir.sampsiz = sampsiz;
o.parameter.dir.strength = strength;
o.parameter.dir.trialnum = trialnum;
o.parameter.dir.sampsiz  = sampsiz;
o.parameter.dir.series = series;
o.parameter.dir.count  = count;
o.parameter.dir.modes =  modes;
o.parameter.dir.mean =  CirCmean.deg.mean;

