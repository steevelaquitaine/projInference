


%SLgetNumPrecision.m
%
%       $Id: SLgetNumPrecision.m $
%        by: steeve laquitaine
%      date: 141107
%   purpose: get numerical precision (number of floating point)
%
%
%     usage:
%
%           precision = SLgetNumPrecision(10.1)

function precision = SLgetNumPrecision(x)

%if integer
precision = 1;

%if not integer
if sum(SLisinteger(x)~=1)>0

    sprintf('(SLgetNumPrecision) Not an integer.')
    
    %precision
    for i = 1 : 10^10
        x = x * 10^i;
        if sum(SLisinteger(x)~=1)==0
            
            precision = 1/10^i;
            break
            
        end
    end
end