function a = auroc(x,y)
%
% AUROC Computes area under ROC curve
%
% Steeve: "The Area Under Curve (AUC) is equal to the probability that a classifier
% will rank a randomly chosen positive instance higher than a randomly
% chosen negative one".
% Trapezoid is a non-parametric method

% -----------
% Input:  x,y: real vectors (not necessarily of same size)
% Output: a: area under ROC curve 
%
% -----------
% 11.08.05 gs
  
  % x: class 0
  % y: class 1
    
  nx = length(x);
  ny = length(y);

  z = flipud(unique([x;y]));
  n = length(z);
  fp = zeros(n+2,1);
  tp = zeros(n+2,1);
  
  fp(1) = 0;
  tp(1) = 0;
  
  for i=1:n
    fp(1+i) = nnz(x>=z(i))/nx;
    tp(1+i) = nnz(y>=z(i))/ny;
  end
  
  fp(end) = 1;
  tp(end) = 1;

  a = trapz(fp,tp);
  