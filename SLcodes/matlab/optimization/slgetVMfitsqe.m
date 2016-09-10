%slgetVMfitsqe(x,y,p)
%
%
% author : steeve laquitaine
%   date : 160115
%purpose : get the squared error (objective function for optimization) 
%          of a least-square fit to von Mises function (not sum to 1)
%
%usage:
%       data = [0 0 0 4 6.8 6.8 4.5 0 0 0];
%       u = 150;
%       k = 40;
%       A = 6.8;
%       [sqe,yfit] = slgetVMfitsqe(1:36:360,data,[u k A])
%       plot(1:36:360,data,'.'); hold on; plot(1:36:360,yfit,'r')

%von mises fit squared error


%von mises fit squared error
function [sqe,yfit] = slgetVMfitsqe(x,y,p)

%get fit parameters
u = p(1);
k = p(2);
A = p(3);

%================
%constrain search
%================
%parameters must be > 0
if any([k;u] < 0)
    sqe = inf;
    yfit = nan;
    return
end
%u must be within 0 and 359
if u > 359
    sqe = 10^300;
    yfit = nan;
    fprintf('%s \n','u was outside the range 0 and 359 deg')
    return
end
%convert to radians
x2rad = SLde2r(x,0);
u2rad = SLde2r(u,0);

%prediction
yfit = A*exp(k.*cos(x2rad-u2rad)-k);

%square error of prediction
sqe = sum((SLmakeRow(yfit) - SLmakeRow(y)).^2);

