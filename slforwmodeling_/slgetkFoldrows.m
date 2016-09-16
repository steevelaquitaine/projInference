


%slgetkFoldrows.m
%
% author : steeve laquitaine
%purpose : returns the row indices of K folds that divide data in each 
%          structure of a Nvar structure 
%           
%usage:
%           f_st,f_end] = slgetkFoldrows(instances,Nf)
%
%inputs :
%
%           instances : a Nvar structure where each structure contains a
%                       matrix of Ni observations x Nd dimensions
%                  Nf : Number of folds in which each structure's matrix
%                       is divided
%
%outputs:
%           f_st : indice of the starting row of each fold x variable
%           f_end: indice of the ending row of each fold x variable

function  [f_st,f_end,NinFold,Nivi] = slgetkFoldrows(instances,Nf)

%number of variables
Nvar = length(instances);

%% for each var (structure)
for i = 1 : Nvar

    %number of instances in for this var
    Nivi = size(instances{i},1);
    
    %number of instances in a fold
    NinFold = round(Nivi/Nf);
    
    %folds indices
    f_sttm = [1 : NinFold : Nivi]';   
    f_st(:,i) = f_sttm;
    f_end(:,i) = [f_st(2:end,i)-1;Nivi];    
end

Nivi = NinFold * Nvar;