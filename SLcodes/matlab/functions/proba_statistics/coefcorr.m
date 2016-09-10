
% steeve laquitaine 14062009
% correlation coefficient + significance

%data: two vector to compare in one matrix 

[r,p] = corrcoef(data,'pairwise')  % Compute sample correlation and p-values.
[i,j] = find(p<0.05);  