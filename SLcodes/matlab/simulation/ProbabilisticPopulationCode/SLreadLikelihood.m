%SLreadLikelihood.m

function o = SLreadLikelihood(o,Optiondisplay)
%logLLH(stimuli|response) is
%   Sum(neurons'responses to stim weighted by log(mean response for this
%   stim))
%
%where the mean response for the stim is defined by the tuning curve
%
%note: can we really safely discard the 2nd
%term in the equation as stated by Jazayeri et al.?
%The second term can be discarded when the
%representation is homogeneous. Is it the case?

o.rs = SLmakeColumn(o.rs);

%likelihood (LLH) from Ma & Jazayeri
A = nan(length(o.diSpace),1);
B = sum(log(gamma(o.rs+1)));
logL = nan(length(o.diSpace),1);
LLH_SeungSompo = nan(length(o.diSpace),1);

for di = 1 : length(o.diSpace)
    
    %Tuning response f at stimulus
    posS = o.diSpace(di);
    fdi = o.f(posS,:)';
    
    %1 - Jazayeri et al.,
    A(di) = sum(fdi);                   %sum over neurons
    %logF  = SLmakeColumn(log(fdi));
    logL(di) = o.rs'*log(fdi) - A(di) - B;
    
    %2 - Exact LLH
    %Seung, & Sompolinski, 1993; Jazayeri et al., 2006; Ma et al.2006
    LLH_SeungSompo(di) = prod(exp(-fdi).* fdi.^(o.rs)./gamma(o.rs+1));
end

%proba (should I do that?)
LLH_Jazayeri = exp(logL)/sum(exp(logL));

%warning
if unique(exp(logL))==0
    fprintf('(SLreadLikelihood) Reaching numerical limits. LogL is too negative because you have: \n')
    fprintf('(SLreadLikelihood) - too many differently tuned neurons "nbNeurons" \n')
    fprintf('(SLreadLikelihood) - too strong tuning width \n')
    fprintf('(SLreadLikelihood) To simulate lower tuning widths, reduce "nbNeurons" \n')
end
LLH_SeungSompo = LLH_SeungSompo./sum(LLH_SeungSompo);

%o
o.LLH_SeungSompo  = LLH_SeungSompo;
o.LLH_Jazayeri = LLH_Jazayeri;

%mean LLH
o.meanLLH_Jazayeri = SLmakeMean(o.diSpace,LLH_Jazayeri);
o.meanLLH_SeungSompo  = SLmakeMean(o.diSpace,LLH_SeungSompo);

%std LLH
o.stdLLH_Jazayeri = SLmakeStd(o.diSpace,LLH_Jazayeri);
o.stdLLH_SeungSompo  = SLmakeStd(o.diSpace,LLH_SeungSompo);

%display
if nargin >1
    o = display3(Optiondisplay,o);
end