% hotelling.m
%
%      usage: p = hotelling(phi1,phi2,r1,r2)
%             phi and r are vectors with the angle in radians
%             and the length (if not set, defaults to all 1s).
%             gives p value that the samples are pulled
%             from the same underlying distribution
%         by: justin gardner
%       date: 12/12/00
%    purpose: Hotelling's two-sample test 
%             
%
function p = hotelling(phi1,phi2,r1,r2,makeplot)

% check arguments
if (nargin == 2)
  phi1 = phi1(find(~isnan(phi1)));
  phi2 = phi2(find(~isnan(phi2)));
  r1 = ones(1,length(phi1));
  r2 = ones(1,length(phi2));
  makeplot = 0;
elseif (nargin == 4)
  good1 = find(~isnan(phi1) & ~isnan(r1));
  good2= find(~isnan(phi2) & ~isnan(r2));
  phi1 = phi1(good1);
  phi2 = phi2(good2);
  r1 = r1(good1);
  r2 = r2(good2);
  makeplot = 0;
elseif (nargin ~= 5)
  help hotelling;
  return
end
  
% convert into column vectors
if (size(phi1,1) == 1), phi1 = phi1';,end
if (size(r1,1) == 1), r1 = r1';,end
if (size(phi2,1) == 1), phi2 = phi2';,end
if (size(r2,1) == 1), r2 = r2';,end

% calculate number of points
n1 = length(phi1);
n2 = length(phi2);

% can't perform test if we don't have enough points
if ((n1 < 2) | (n2 < 2))
  p = nan;
  return;
end

% find the x and y's of the points
x1 = r1.*cos(phi1);
y1 = r1.*sin(phi1);
x2 = r2.*cos(phi2);
y2 = r2.*sin(phi2);

% calculate the means
m1 = mean([x1 y1]);
m2 = mean([x2 y2]);

% calculate the covariance matrix
c1 = cov([x1 y1]);
c2 = cov([x2 y2]);

% calculate pooled sum of squares
ssx = c1(1,1)*(n1-1) + c2(1,1)*(n2-1);
ssy = c1(2,2)*(n1-1) + c2(2,2)*(n2-1);

% caluclate correlation coefficent
r = (c1(1,2)*(n1-1)+c2(1,2)*(n2-1))*(ssx*ssy)^(-1/2);

% calculate t-statistic
t1 = (m1(1)-m2(1))*(((1/n1)+(1/n2))*ssx/(n1+n2-2))^(-1/2);
t2 = (m1(2)-m2(2))*(((1/n1)+(1/n2))*ssy/(n1+n2-2))^(-1/2);

% calculate T squared value
T2 = ((1-r^2)^-1)*(t1^2-2*r*t1*t2+t2^2);

% convert T2 to an F value
F = T2*(n1+n2-3)/(2*(n1+n2-2));

% calculate p-value from F distribution
% using betainc function (see. Numerical Recipies p.181)
v1 = 2;v2 = n1+n2-3;
p = betainc(v2/(v2+v1*F),v2/2,v1/2);

% plot the data
if (makeplot)
  figure
  % draw data and 95% confidence ellipses
  rconf(phi1,r1,.95,'g');
  rconf(phi2,r2,.95,'b');
  axis([-max(abs([x1 ; x2])) max(abs([x1 ; x2])) ...
	-max(abs([y1 ; y2])) max(abs([y1; y2]))]);
  zoom on
end








