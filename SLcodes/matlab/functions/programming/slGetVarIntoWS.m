

%slGetVarIntoWS.m
%author: steeve laquitaine
%purpose: make workspace variables accessible to function
%
%usage: just put into function 
%
%          slGetVarIntoWS.m
%
%ref
%   http://www.mathworks.com/matlabcentral/answers/719-declaring-many-variables-from-work-space-into-function-workspace


%Get all base workspace variables into this function.
T = evalin('base','whos');
for ii = 1:length(T)
    C_ =  evalin('base',[T(ii).name ';']);
    eval([T(ii).name,'=C_;']);
end