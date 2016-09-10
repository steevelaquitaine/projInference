


%author: steeve laquitaine
%purpose : check whether a number y is located between two numbers x1 and
%          x2

function o = slIsBetween(y,x1,x2)

o = 0;
if y >= x1 && y <= x2
    o = 1;
end
    
    