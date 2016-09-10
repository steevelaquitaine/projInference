

%y and x are a column vector


function CIByCond = slCalculateCISortedBy1Var(data,x)

%locate each condition
[conds,~,cond_location] = unique(x);

%average data for each condition
meangByCond = accumarray(cond_location,data,[length(conds) 1],@mean);
emargByCond = accumarray(cond_location,data,[length(conds) 1],@makeCI95);
CIByCond = [meangByCond-emargByCond meangByCond+emargByCond];
             
%calculaye error margin for confidence interval
function err_margin = makeCI95(data)

m = nanmean(data);
s = nanstd(data);
n = length(data);
err_margin = 1.96*s/sqrt(n);
