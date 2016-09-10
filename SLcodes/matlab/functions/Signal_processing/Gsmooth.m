function [ output_args ] = Gsmooth(h, w)
%GSMOOTH Summary of this function goes here
%   Detailed explanation goes here

% Algorithm from neuroexplorer (steeve laquitaine 27/08/2011)
% The parameters of the filter are such that the width of the
% Gaussian curve at half the height is the specified filter width.
% Please note that for the Gaussian filter, width of the filter (w)
% can be non-integer.  For example, w can be 3.5.

% Input
    %- h: time serie to smooth
    % -w: window.
    
% Output:
    % sh: smoothed time serie with a gaussian filter
    
% The Gaussian filter is calculated using the following formula:
d = ( w + 1 )/2;
sigma = - w * w * 0.25 / log(0.5);
i=-2*d : 2*d;
norm = sum(exp( -i .* i / sigma));
f= exp(-i.*i/sigma)/norm;
f=f';

h2=nan(numel(i), numel(h));
for j=2*d+1 :numel(h)-2*d;
    h2(1:numel(i),j)=h(j-2*d: j +2*d);
end

% repeat the datas at the extremity
h2(:,1:2*d)=h2(:,2*d+1:4*d);
h2(:,numel(h)-2*d+1:end)=h2(:,   numel(h)+1 -4*d:numel(h)-2*d);

% and the smoothed histogram sh[i] is calculated as   sh[i] =
% ....f[j-1]*h[i-1] + f[0]*h[i] + f[1]*h[i+1]...

for k=1: numel(h)
   sh(k)=sum(h2(:,k).* f);
end

end

