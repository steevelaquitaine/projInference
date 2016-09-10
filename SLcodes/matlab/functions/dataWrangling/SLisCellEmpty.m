

% author: steeve laquitaine
%   date: 15/04/01
%purpose: get empty cells

function emptyCells = SLisCellEmpty(data)

emptyCells = cellfun('isempty',data);