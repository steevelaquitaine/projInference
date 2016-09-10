

%SLaddPoissonNoise.m
%
% author: steeve laquitaine
%   date: 150631
%purpose: add poisson noise the
%
%usage: 
%
%   o.nNeuron = 40;
%   o.rSpace = 0 : 1 : 60;
%   o.diSpace = 1 : 1 : 360;
%   o.TuningDiSpace = o.diSpace;
%   o.f % a D directions * N neurons matrix of mean neural responses to
%       directions
%
%   o = SLaddPoissonNoise(o,Optiondisplay)

%Poisson noise 
function o = SLaddPoissonNoise(o,Optiondisplay)

%check
if ~isfield(o,'rSpace')
    fprintf('(SLaddPoissonNoise) I am missing o.rSpace the range of firing responses (e.g., o.rSpace = [0:1:60]) \n')
    keyboard
end
if ~isfield(o,'nNeuron')
    fprintf('(SLaddPoissonNoise) I am missing o.nNeuron the number of neurons (e.g., o.nbNeuron = 360) \n')
    keyboard
end
if ~isfield(o,'f')
    fprintf('(SLaddPoissonNoise) I am missing o.f the neurons tuning vectors \n')
    keyboard
end

%check that tuning and LLH are defined on same space
if sum(o.diSpace - o.TuningDiSpace) ~ 0;
    fprintf('(SLaddPoissonNoise) "o.diSpace" must = "o.TuningDiSpace" \n')
end

%P(responses|direction Di ,neuron Ni) is defined by a Poisson
%distribution.
%note: for non integer values the gamma function interpolates
%factorial function. Factorial(n) = gamma(n+1);
o.PRgivenDi = nan(numel(o.rSpace),length(o.diSpace),o.nNeuron);

%vectorize rSpace for speed
rSpaceAll = o.rSpace(ones(length(o.diSpace),1),:,ones(1,o.nNeuron));
rSpaceAll = permute(rSpaceAll,[2 1 3]);

%vectorize o.f for speed
fAll = o.f(:,:,ones(1,length(o.rSpace)));
fAll = permute(fAll,[3 1 2]);

%P(responses|direction Di ,neuron Ni) is defined by a Poisson
%distribution.
%note: for non integer values the gamma function interpolates
%factorial function. Factorial(n) = gamma(n+1);
o.PRgivenDi = (exp(-fAll) .* fAll.^rSpaceAll)./gamma(rSpaceAll+1);
o.PRgivenDi = bsxfun(@rdivide,o.PRgivenDi,sum(o.PRgivenDi,1)); %proba

%check if display
if nargin > 1
    o = display1(Optiondisplay,o);
end
