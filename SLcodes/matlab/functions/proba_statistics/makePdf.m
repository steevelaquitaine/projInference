
%makePdf.m
%
% author: steeve laquitaine
%   date: 140805
%purpose: make a probability distribution
%
%  usage:
%        [p,xpdf,c] = makePdf(5:10:360,10*randn(10,1)+180,'raw')
%        [p,xpdf,c] = makePdf(5:10:360,motdir,'smooth')
%

function [p,xpdf,c] = makePdf(x,dataset,type)

c = [];

%case raw pdf
if strcmp(type,'raw')
    
    %create pdf (raw)
    c = histc(dataset,x);
    
    %probability
    %warning
    if c==0
        fprintf('(makePdf) p=0 ')
        keyboard
    end
    p = c/sum(c);
    xpdf = x;
end

%case smooth pdf
if strcmp(type,'smooth')
    xpdf = min(x):1:max(x);
    p = ksdensity(dataset,xpdf);
    
    %warning
    if p==0
        fprintf('(makePdf) p=0 ')
        keyboard
    end
    p = p/sum(p);
    

end
