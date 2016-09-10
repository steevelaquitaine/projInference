function [Power] = powerAOVI(F,v1,v2,m,a)
%  Statistical Power of a Performed ANOVA Test Model I 
%
%   Syntax: function [Power] = powerAOVI(F,v1,v2,m,a) 
%      
%     Inputs:
%          F - observed F-statistic value. 
%         v1 - numerator degrees of freedom (k - 1; k = number of treatments). 
%         v2 - denominator degrees of freedom (N - k; N = total number of data).
%          m - method of noncentrality parameter estimate (1=unbiased;
%              else=uniformly minimum variance unbiased)
%          a - significance level. 
%     Output:
%          Power - the output power of the performed ANOVA test Model I.
%          Graphic - Central and noncentral F-distributions of the ANOVA Model I
%  
%    Example: For a balanced design with k = 3, n = 10 (N = k*n), F = 4.60375, then
%             v1 = 2 and v2 = 27. We are interested to get the post-hoc power with
%             a nonecetrality parameter uniformly minimum variance unbiased estimate,
%             and a significance level = 0.05.
%
%             Calling on Matlab the function: 
%                powerAOVI(4.60375,2,27,2,0.05)
%             Answer is:
%                ans= 0.5645
%           
%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  Copyright (C) November 2002.
%        $Revised July 2008.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2002). powerAOVI: Calculates the 
%    statistical power of a performed ANOVA test Model I. A MATLAB file.
%    [WWW document]. URL http://www.mathworks.com/matlabcentral/fileexchange/
%    loadFile.do?objectId=2696
%
%  References:
% 
%  Burden, R. and Douglas, F. (1985), Numerical Analysis. PWS Boston.
%  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995). Continuous Univariate
%              Distributions, Vol. 2, 2nd edition. New York: John Wiley.
%  Winer, B.J. (1971), Statistical Concepts in Experimental Design.
%              New York: McGraw-Hill. pp. 220-222; 225-228.
%  Zar, J.H. (1999), Biostatistical Analysis (2nd ed.).
%              NJ: Prentice-Hall, Englewood Cliffs. p. 193. 
% 

if nargin < 5, 
    a = 0.05; 
end 

if ~isscalar(a)
   error('POWERAOVI requires a scalar ALPHA value.');
end

if ~isnumeric(a) || isnan(a) || (a <= 0) || (a >= 1)
   error('POWERAOVI requires 0 < ALPHA < 1.');
end

if nargin < 4, 
    m = 2; 
end 

if nargin < 3, 
    error('Requires at least three input arguments.'); 
end 

p=1-a;
sl=a;
%
%The deeply framework of the analysis of variance fundamentals for Model I (fixed effects),
%establish a very good knowledgement of the factor to be tested (independent variable) or
%controled part of the total observed variability, and one would expecting the variability
%due to that factor should be greater than the variability of the random error or not 
%controled part. So, the observed F-statistic value should be greater than 1.0 and a
%noncentrality parameter greater than zero. Otherwise, the noncentrality parameter results
%to be a negative or inexistent value. However, this Matlab file considering a such a case
%displaying a WARNING and consequently the statistical power value equals the significance
%level.
%
if (F < 1.0)
   disp('WARNING: For a true noncentrality value the observed F-statistic value should be greater than 1.0');
   l = 0.0;
else
    if (m == 1),
        disp('Based on the unbiased noncentrality parameter estimate.');
        l = v1*(F-1); %Estimated noncentrality parameter from the observed F-statistic value.
    else
        disp('Based on the uniformly minimum variance unbiased noncentrality parameter estimate.');      
        l = v1*(F-1);
        l = ((l+v1)*(v2-2)/v2)-v1; %Based on the uniformly minimum variance unbiased estimate
                                 %of non-centrality parameter (Johnson et al., 1995)
    end
end
%
% Aproximation of the noncentral F distribution to central F distribution.
%
finv(p,v1,v2);
Fc = finv(p,v1,v2)/(1+(l/v1)); %Expected F-statistic value adjusted to the estimated
%noncentrality parameter.
v1m = ((v1+l)^2)/(v1+(2*l)); %Numerator degrees of freedom adjusted to the estimated
%noncentrality parameter.
%
% Because the numerator degrees of freedom corrected by the noncentrality parameter
% could results a fraction, the probability function associated to the F distribution
% function is resolved by the Simpson's 1/3 numerical integration method.
%
x = linspace(.00001,Fc,10001);
DF = x(2)-x(1);
y = ((v1m/v2)^(.5*v1m)/(beta((.5*v1m),(.5*v2))));
y = y*(x.^((.5*v1m)-1)).*(((x.*(v1m/v2))+1).^(-.5*(v1m+v2)));
N1 = length(x);
Power = 1-(DF.*(y(1)+y(N1) + 4*sum(y(2:2:N1-1))+2*sum(y(3:2:N1-2)))/3.0);
F = finv(1-a,v1,v2);
x = (0.01:0.1:10.01)';x1=F;
p1 = ncfpdf(x,v1,v2,l);
p = fpdf(x,v1,v2);y=(0.0:.001:max(p)+.1);
plot(x,p1,'--k')
hold on
plot(x,p,'-k')
D = [x p1];
y1 = interp1(D(:,1),D(:,2),x1);
g = find(D(:,1)>x1);
X = D(:,1);
Y = D(:,2);
X = [x1;X(g)];
Y = [y1;Y(g)];
hold on
a = area(X,Y);
set(gcf,'Renderer','OpenGL')
set(a,'FaceColor',[0.75,0.75,0.75])
hold on
x2 = [x1,x1];
y2 = [y1,.5];
plot(x2,y2,'-.k')
xlabel('\itF\rm-ratio');
ylabel('Frequency of \itF\rm-ratios');
text(F+.05,0.4,['\leftarrow \itF_{',num2str(sl) ', ',num2str(v1m) ', ',num2str(v2) '}\rm = ',...
    num2str(F)],'FontWeight','bold','FontSize',9)
hold on
text(F+F+.32,0.52,['Power = ',num2str(Power)],'FontWeight','bold')
hold on
rectangle('Position',[F+F,.5,.24,0.03],'Curvature',[0,1],'FaceColor',[0.75,0.75,0.75])
if (m == 1),
    title({'Central (-) and noncentral (- -) \itF\rm-distributions of the ANOVA Model I based on';...
        'the unbiased noncentrality parameter estimate.'}),
else
    title({'Central (-) and noncentral (--) \itF\rm-distributions of the ANOVA Model I based on';...
        'the uniformly minimum variance unbiased noncentrality parameter estimate.'}),
end
hold off

return,