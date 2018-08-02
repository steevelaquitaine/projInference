% SLcircConv.m
%
%
%   date: 24/11/2013
% Author: Steeve Laquitaine
%purpose: Calculate circular convolutions with circular data
%         The probability that value i in vector 2 would be combined with at least one value from vector 1.
%
%  Usage:
%
%     cconv = circConv(vector 1,vector 2)
%     vector 1 and 2 are col vectors (vertical) or matrices for convolution 
%     column by column. e.g., two probability densities.
%     
%references: http://www.mathworks.com/help/signal/ug/linear-and-circular-convolution.html
%
%note: convolutions are distributive
%e.g., case1=case2
%case1 = circConv(0.1,0.2)+circConv(0.3,0.2)
%case2 = circConv(0.1+0.3,0.2)

function cconv = SLcircConv(v1,v2)

cconv = ifft(fft(v1).*fft(v2));

%sanity check
if any(~isreal(cconv))
    fprintf('%s \n','Convolution failed because of not real values (e.g., NaN) in v1 or v2')
    dbstack
    keyboard
end
    

