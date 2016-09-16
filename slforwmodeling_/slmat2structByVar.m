
%slmat2structByVar.m
%
%
%author : steeve laquitaine
%purpose: sort matrix data where matrix rows are data associated with 
%         different repeated data in a N variables structure 
%
%usage :
%
%   [M_struc,v_u] = slmat2structByVar(rand(100,5),repmat([1:5]',20,1))
%
%inputs
%
%       M : Ni x Nv matrix where each row is associated with a variable in
%           v
%       v : Ni vector of variable
%
%outputs: 
%
% M_struc : M sorted in a structure for each unique variable value
%     v_u : unique variable values
%
function [M_struc,v_u] = slmat2structByVar(M,v)

%unique variables
v_u = unique(v);

%number of unique variables
Nvi = length(v_u);

for i = 1 : Nvi
    M_struc{i} = M(v==v_u(i),:);    
end
