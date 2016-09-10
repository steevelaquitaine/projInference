

%slmakeCircCI.m
%
% author: steeve laquitaine
%purpose: calculate confidence intervals over circular data
%
%  usage: 
%
%          ci = slmakeCircCI(SLde2r([0 90 180 90],1),0.95)
%
%
%reference : Statistical analysis of circular data, N. I. Fisher


function ci = slmakeCircCI(angles,level)

angles = SLmakeColumn(angles);

%variable for confidence intervals
%calculate mean resultant vector length
n = length(angles);
r = abs(sum(exp(1i*angles)))/n;
R = n.*r;
c2 = chi2inv(level,1);

% check for resultant vector length and select appropriate formula
ci = zeros(size(r));

for i = 1 : numel(r)
  if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
    ci(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));
  elseif r(i) >= .9
    ci(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      
  end
end

ci = acos(ci./R);


