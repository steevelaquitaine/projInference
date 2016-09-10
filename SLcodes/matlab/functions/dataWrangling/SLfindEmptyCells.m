

%SLfindEmptyCells.m
%
%       $Id: SLgetAxes.m $
%        by: steeve laquitaine
%      date: 140709
%   purpose: get figures axes
%
%     usage:
%
%           emptyCells = SLfindEmptyCells([{1},{''},{2}])
%

function emptyCells = SLfindEmptyCells(x)
emptyCells = cellfun(@isempty,x);
