function p=nrm_colsmooth_ntimes(x)
% steeve last version 12/ 10/2008 

% use a matrix
% smooth the matrix columns with a gaussian moving window normalised over
% the vector length n times

clc;

% Enter parameter
n=100;

n=n-1;
p=m_nanngmaw_col(x);
% n times
for i=1:n
    p=m_nanngmaw_col(p);
end


