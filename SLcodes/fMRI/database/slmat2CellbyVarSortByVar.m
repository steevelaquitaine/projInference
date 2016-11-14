
%slmat2CellbyVarSortByVar.m
%
%
%author : steeve laquitaine
%purpose: sort matrix data where matrix rows are data associated with 
%         different repeated data in a N variables structure 
%
%  usage:
%
%   [M_struc,v_u,M_ix] = slmat2CellbyVarSortByVar(rand(100,5),repmat([1:5]',20,1),repmat([1:2]',50,1),1)
%
%inputs
%
%       M : Ni x Nv matrix where each row is associated with a variable in
%           v
%      v1 : Ni vector of variable 1
%  v2_all : Ni vector of variable 2
%  v2_in  : value of variable 2 used to sub-select data
%
%
%outputs: 
%
% M_struc : M sorted in a structure for each unique variable value
%     v_u : unique variable values
%
function [M_struc,v1_u,M_ix] = slmat2CellbyVarSortByVar(M,v1,v2_all,v2_in)

%unique variables
v1_u = unique(v1);

%number of unique variables
Nvi = length(v1_u);

for i = 1 : Nvi
    M_ix{i} = find(v1==v1_u(i) & v2_all == v2_in);
    M_struc{i} = M(M_ix{i},:);    
end
