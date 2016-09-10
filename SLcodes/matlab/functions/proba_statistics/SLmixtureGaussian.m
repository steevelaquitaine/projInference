

%SLmixtureGaussian.m
%
% author: Steeve Laquitaine
%   date: 150109
%purpose: Create a mixture of Gaussian probability distribution
%
%
%usage:
%
%       MoG = SLmixtureGaussian(5:10:355,[145 305],[20 20])
  

function MoG = SLmixtureGaussian(x,modes,stds)

f = 0.5*(1./(stds(1)*sqrt(2*pi))) * exp( -.5 * ((x - modes(1))./stds(1)).^ 2 ) + ...
    0.5*(1./(stds(2)*sqrt(2*pi))) * exp( -.5 * ((x - modes(2))./stds(2)).^ 2 );

%probability
MoG = f./sum(f);
