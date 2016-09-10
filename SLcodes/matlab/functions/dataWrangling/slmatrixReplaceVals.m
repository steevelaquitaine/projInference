

%author : steeve laquitaine
%  date : 150726
%purpose : replace values in a matrix by new values
%
%   usage: 
%
%       xnew = slmatrixReplaceVals([1 2; 1 3],[4 5 6])

function xnew = slmatrixReplaceVals(xold,rep)

xoldvec = unique(xold,'stable');
nxoldvec = length(xoldvec);

xnew = nan(size(xold));

for i = 1 : nxoldvec 
    
    xnew(xold==xoldvec(i)) = rep(i);
    
end