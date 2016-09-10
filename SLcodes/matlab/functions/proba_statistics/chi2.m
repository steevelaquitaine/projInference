function [chi,p]=chi2(tab);

% [chi,p]=chi2(tab);
% 
% Chi-square test for a n rows x 2 columns contingency table.
%
% Input
% tab : matrix of integer (n rows and 2 columns) giving simple frequencies
% (typically = number of occurrence of a binomial variable like
% "rain/no rain" across a family of index (like weather types)).
%%	       Col 1     Col 2	   
%Ligne 1	3	      20	   
%Ligne 2	12	      12	   
%Ligne n	19	      15	
% Output
% 'chi' ; scalar giving the value of chi-square
% 'p' : scalar giving the level of significance
%
% Vincent Moron
% Nov. 2005

[row,column]=size(tab);
t1=sum(tab);
t2=sum(tab');
t=sum(tab(:));
for i=1:row;
    e(i,1)=(t1(1)*t2(i))./t;
    e(i,2)=(t1(2)*t2(i))./t;
end

tab=tab(:);
e=e(:);
chi=sum(((tab-e).^2)./e);
dof=row-1;
p=1-chi2cdf(chi,dof); 



