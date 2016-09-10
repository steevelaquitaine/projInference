
% SLfindRow.m
%
%     author: steeve laquitaine
%       date: 140515 last modif 140808
%    purpose: find the index of a row in a matrix
%      usage:
%         indx = SLfindRow([1 2],[5 6; 1 2;3 4; 1 2]) 

function indx = SLfindRow(rowValues,Matrix)

%case rowValue is larger than Matrix. Not possible.
if size(rowValues,1) > size(Matrix,1)
    fprintf('%s \n', 'It s the other around. Matrix is 1st arg then rowValue.')
    return
end

[~,indx] = ismember(Matrix,rowValues,'rows');
