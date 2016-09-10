

%author: steeve laquitaine
%usage :
%
%       p = SLdeltaPdf(1:1:360,[90 95],0.24,1)

function [p,x] = SLdeltaPdf(x,mode,pMode,maxProba)
 
%initialize
nX = length(x);
nMode = length(mode);
p = nan(nX,nMode);

%probability non mode
p(1:nX,:) = repmat((maxProba - pMode)./nX,nX,nMode);
p(bsxfun(@ne,x',mode)) = (maxProba - pMode)./nX;

%warning
if isempty(intersect(mode,x))
    fprintf('(SLdeltaPdf) The mode must be included in x')
end
    
%probability of mode
p(bsxfun(@eq,x',mode)) = pMode;

%warning
if sum(p)<0.990
    fprintf('(SLdeltaPdf) sum(p) =~ 1. Please check the result \n')
end
    