
% SLKtoStdev.m
%
%     Author: Steeve Laquitaine
%    purpose: calculate standard deviation of a von Mises
%       
%       usage: 
%
%           S = SLKtoStdev(25);
%               
%  adjusted for large k.

function S = SLKtoStdev(k)

%setup
x = 1:1:360;
u = 180;
p = vmPdfs(x,u,k,'norm');

%angular distance
dis = SLmakeColumn(SLvectors2signedAngle(x,u,'polar'));
dis = dis(:,ones(numel(k),1));

%std
S = sqrt(sum(p.*(dis.^2)));