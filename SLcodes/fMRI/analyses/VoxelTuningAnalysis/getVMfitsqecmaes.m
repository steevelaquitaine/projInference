
%von mises fit squared error
function [sqe,yfit] = getVMfitsqecmaes(p,x,y)

%----------------- vMises ------------------------
%get fit parameters
u=p(1);
k=p(2);
A=p(3);

%parameters must be > 0
if any([k;u] < 0)
    sqe = inf;
    yfit = nan;
    return
end
    
%convert to radians
x2rad = SLde2r(x,0);
u2rad = SLde2r(u,0);

%prediction
yfit = A*exp(k.*cos(x2rad-u2rad)-k)./(2*pi.*besseli(0,k,1));

%square error of prediction
sqe = sum((yfit' - y).^2);
