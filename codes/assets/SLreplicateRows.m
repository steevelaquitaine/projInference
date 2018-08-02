
% SLreplicateRows.m
%
%     author: steeve laquitaine
%       date: 140620
%    purpose: repeat rows of a matrix
%
%      usage:
%           RM=SLreplicateRows([0 0; 2 2; 3 3],3)
%
%references
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/239627

function RM=SLreplicateRows(M,Rep)

%call for help (may slow down code)
% if ieNotDefined('M')
%     help SLreplicateRows
%     return
% end

nRows=size(M,1);
RepAllRows=Rep(1,ones(1,nRows));
idx = accumarray((cumsum(RepAllRows)+1)',ones(size(RepAllRows))');
RM = M(cumsum(idx(1:end-1))+1,:);