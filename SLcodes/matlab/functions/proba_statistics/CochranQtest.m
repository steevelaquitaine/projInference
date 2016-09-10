function CochranQtest(X,alpha)
%COCHRANQTEST Cochran's Q Test for Margins Homogeneity.
%   COCHRANQTEST performs the Cochran's Q test for multi-way tables which
%   each variable has two levels, that is, for 2x2...x2 tables. It is used
%   to test the homogeneity of the one-dimensional margins. So, Cochran's 
%   Q test is a nonparametric test examining change in a dichotomous variable
%   across more than two observations. It can be approximate to a Chi-square 
%   statistic. When there are only two binary response variables, Cochran's 
%   Q test simplifies to McNemar's test. The test was proposed by Cochran (1950).
%
%   Then, for a NxK table where N is the number of subjects or blocks and K
%   is the number of repeated measures or different tratments. The null 
%   hypothesis to test is,
%
%         Ho: p_1 = p_2 = . . . = p_K; j = 1,2,...,K.
%
%   p_j is the probability of a success for a case under condition/treatment j.
%
%   Syntax: function CochranQtest(X,alpha) 
%      
%   Input:
%         X - data matrix (Size of matrix must be N-by-K;binary value=column 1,
%             treatment or repeated measure=column 2;subject or block=column 3). 
%     alpha - significance level (default = 0.05).
%
%   Output:
%         A table with:
%         - Cochran statistic, number of subjects or blocks, number of repeated
%           mesures, degrees of freedom and upper-tail P-value.
%
%   Example: From the hypothetical example given in the internet site URL 
%            http://www.bof.fire.ca.gov/pdfs/Lewis-HMP.pdf on Section 3-7. If we 
%            observe a binary response (say, effective (=1) and not effective (=0) 
%            at 8 locations (N) on 4 occasions (K), we might be interested in testing
%            with a significance level = 0.05 if the effectiveness rate changes over
%            time. 
%            
%                                             Time
%                   ----------------------------------------------
%                    Location       1       2       3       4  
%                   ----------------------------------------------
%                        1          1       1       0       1
%                        2          1       1       1       1            
%                        3          0       0       0       0
%                        4          0       0       0       0
%                        5          1       1       1       1
%                        6          1       0       0       0
%                        7          1       1       1       1        
%                        8          1       1       1       0       
%                   ----------------------------------------------
%                                       
%   Data matrix must be:
%      X=[1 1 1;1 1 2;0 1 3;0 1 4;1 1 5;1 1 6;1 1 7;1 1 8;
%      1 2 1;1 2 2;0 2 3;0 2 4;1 2 5;0 2 6;1 2 7;1 2 8;
%      0 3 1;1 3 2;0 3 3;0 3 4;1 3 5;0 3 6;1 3 7;1 3 8;
%      1 4 1;1 4 2;0 4 3;0 4 4;1 4 5;0 4 6;1 4 7;0 4 8];
%
%   Calling on Matlab the function: 
%             cochranqtest(X)
%
%   Answer is:
%
%   Table for the Cochran's Q test.
%   -------------------------------------------------------
%     Q statistic       N         K        df         P  
%   -------------------------------------------------------
%        3.6667         8         4         3       0.2998
%   -------------------------------------------------------
%   If the P-value is smaller than 0.05
%   the Ho tested results statistically significant. Otherwise,
%   it is not significative.
%  
%  Created by A. Trujillo-Ortiz, R. Hernandez-Walls and A. Castro-Perez
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%   Copyright. November 12, 2004.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A., R. Hernandez-Walls and A. Castro-Perez. (2004). CochranQtest:
%    Cochran's Q Test for Margins Homogeneity. A MATLAB file. [WWW document]. URL 
%    http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6431
%
%   References:
%   Cochran, W. G. (1950), The comparision of percentage in matched 
%            samples. Biometrika, 37:256-266.
%

if nargin < 2,
   alpha = 0.05;  %(default) 
elseif (length(alpha)>1),
   error('Requires a scalar alpha value.');
elseif ((alpha <= 0) | (alpha >= 1)),
   error('Requires 0 < alpha < 1.');
end;

if nargin < 1, 
   error('Requires at least one input argument.');
   return;
end;

if  min(X(:,1)) ~= 0 | max(X(:,1)) ~= 1
    error('Requires an input binary variable: 0 or 1.');
end;

k = max(X(:,2));
s = max(X(:,3));

M = sum(X(:,1))/length(X(:,1));

Mt=[];
indice = X(:,2);
for i = 1:k
   Xe = find(indice==i);
   eval(['X' num2str(i) '=X(Xe,1);']);
   eval(['mt' num2str(i) '=mean(X' num2str(i) ');']);
   eval(['x = mt' num2str(i) ';']);
   Mt = [Mt;x];
end;

T = [];
for i = 1:k
    eval(['x = (Mt(' num2str(i) ') - M).^2;']);
    T = [T;x];
end;

Mn=[];
indice = X(:,3);
for j = 1:s
   Xe = find(indice==j);
   eval(['X' num2str(j) '=X(Xe,1);']);
   eval(['mn' num2str(j) '=mean(X' num2str(j) ');']);
   eval(['x = mn' num2str(j) ';']);
   Mn = [Mn;x];
end;

N = [];
for j = 1:s
    eval(['x = (Mn(' num2str(j) ')*(1 - Mn(' num2str(j) ')));']);
    N = [N;x];
end;

Num = sum(T);
Den = sum(N);

Q = (s^2*(k - 1)/k)*Num/Den;
v = k - 1;

P = 1 - chi2cdf(Q,v);    

disp(' ')
disp('Table for the Cochran''s Q test.')
fprintf('-------------------------------------------------------\n');
disp('  Q statistic       N         K        df         P  '); 
fprintf('-------------------------------------------------------\n');
fprintf(' %10.4f%10.i%10.i%10.i   %10.4f\n',[Q,s,k,v,P].');
fprintf('-------------------------------------------------------\n');
fprintf('If the P-value is smaller than% 3.2f\n', alpha );
disp('the Ho tested results statistically significant. Otherwise,');
disp('it is not significative.');

return;