function [y,sorted_x ]=find_betw(x,mini, maxi,idx)
    
% steeve laquitaine 
% 08/06/2011 : 18:07

% DETAILED EXPLANATION GOES THERE

if num2str(idx)=='strict';
    y=find(x < maxi & x > mini);
else
    y=find(x<=maxi & x >= mini);
end

sorted_x=x(y);

end
