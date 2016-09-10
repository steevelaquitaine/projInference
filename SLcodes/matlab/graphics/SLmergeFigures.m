
%SLmergeFigures.m
%
%
% author: steeve laquitaine
%purpose: merge plotted figures
%
%  usage:
%       SLmergeFigures([3 4])


function SLmergeFigures(myfigures)

%get figures opened
hf = SLgetFigures;
hf = sort(hf);
hf = hf(myfigures);

%copy together on a another figure
hf(end+1) = figure('color','w');
for i = 1 : length(hf) - 1
    
    %get all plot objects
    hc  = get(hf(i),'children');
    hgc = get(hc, 'children');
    
    %make cell
    if ~iscell(hgc)
        hgc = mat2cell(hgc);
    end
    
    %copy paste on new figure
    for j = 1 : length(hgc)
        copyobj(hgc{j},gca);
    end
end
