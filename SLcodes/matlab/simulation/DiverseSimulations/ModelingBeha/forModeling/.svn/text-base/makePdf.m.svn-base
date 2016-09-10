% Calculate the smoothed pdf of a dataset
function [p xpdf]=makePdf(x,dataset,type)

if strcmp(type,'raw')
    % create pdf (raw)
    c=histc(dataset,x);
    p=c/sum(c);
    xpdf=x;
end

if strcmp(type,'smooth')
    % create pdf (smoothed)
    xpdf=min(x):1:max(x);
    p=ksdensity(dataset,xpdf);
    p=p/sum(p);
end
