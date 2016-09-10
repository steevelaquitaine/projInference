

% author: steeve laquitaien
%purpose: add zero before a number when converting to string
%
%usage
%
%       x = sladdZerobeforeNum2str(1)

function x = sladdZerobeforeNum2str(x)

if x < 10
    x = ['0' num2str(x)];
else 
    x = num2str(x);
end
