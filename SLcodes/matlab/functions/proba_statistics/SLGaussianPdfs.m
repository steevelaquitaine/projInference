
%SLgaussianPdfs.m
%
%      author: steeve laquitaine
%        date: 141105
%     purpose: compute a 1 dimensional (x) Gaussian
%
%       usage: 
%
%              x = -100:1:100; 
%              p = SLGaussianPdfs(x,0,1)
%
%
%Description:
%
%       x = space
%       u = mean of x
%       stdev = stdev of x
%       type: 'norm' if you want probability that
%              sum to 1.
%
%
function p = SLGaussianPdfs(x,u,stdev)

%Gaussian (sum to 1)
dist = x - u;
scale = 1./(sqrt(2*pi)*stdev);
p = scale.*exp(-0.5*(dist./stdev).^2);
p = p';
