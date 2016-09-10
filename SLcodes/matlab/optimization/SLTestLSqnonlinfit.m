
%SLTestLSqnonlinfit.m
%
%author: steeve laquitaine
%  date: 140804
%purpose: example of usage of Least square non linear fit solver. Find the
%best parameter a that help predict y,
%
%  usage:
%       [x,resnorm,residual,exitflag,output,lambda,jacobian] = SLTestLSqnonlinfit;
%

function [x,resnorm,residual,exitflag,output,lambda,jacobian] = SLTestLSqnonlinfit;

%Variable to predict
y = SLmakeColumn(1:10);

%Model's initial parameters (guesses)
a0 = 2;

%Fit
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@(a) makeError(y,a), a0);

%Error (residual)
function error = makeError(y,a)

%variable
x = SLmakeColumn(1:10);

%model
pred = a.*x + 0;

%error
error = y - pred;