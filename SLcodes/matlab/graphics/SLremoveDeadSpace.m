%SLremoveDeadSpace.m
%
%       $Id: SLremoveDeadSpace.m $
%        by: steeve laquitaine
%      date: 140528 last modif. 140806
%   purpose: remove dead space around the axes of figure
%
%   usage: 
%             figure;
%             count=0;
%             for i=1:1
%                 for j=1:2
%                     count=count+1;
%                     subplot(1,2,count);
%                 end
%             end
%             SLremoveDeadSpace
%
%             SLremoveDeadSpace(0.05)
%
%need to solve that sometimes cut figure which is bad !

function SLremoveDeadSpace(varargin)

%default percDeadSpace that always works nice.
%percDeadSpace = 0.12; %old value
if ~isempty(varargin)
    percDeadSpace = varargin{:};
else
    percDeadSpace = 0.05;
end

%figure's axes
allAxesInFigure = SLgetAxes;

%remove legends handles
%----------------------
lgH = SLgetLegendAxes;
allAxesInFigure = setdiff(allAxesInFigure,lgH);
[allAxesInFigure,~] = sort(allAxesInFigure,'ascend');

%check XLabel
%------------
%if Xlabel, add deadSpace between rows
numAxes = length(allAxesInFigure);

%case nb axes > 1
if length(allAxesInFigure) > 1
    allXLbs = cell2mat(get(allAxesInFigure,'xlabel'));
    XLbPos = cell2mat(get(allXLbs,'extent'));
    pos = cell2mat(get(allAxesInFigure,'position'));
else
    %case 1 axis
    allXLbs = get(allAxesInFigure,'xlabel');
    XLbPos = get(allXLbs,'extent');
    pos = get(allAxesInFigure,'position');
end
nCols = numel(unique(pos(:,1)));
nRows = numel(unique(pos(:,2))); 

%check if there are xlabels
%count how many rows have XLabels
isXlabel = SLisXLabel(allAxesInFigure);
[~,~,rowsID] = unique(pos(:,2),'stable');

%check if Xlabel exists for each row. If yes get max XLabel height 
%for this row and Xlabels all axes of this same row to place the axes at
%the same level.
for i = 1 : nRows
    
    %check if Xlabel exists for this row
    thisRow = i==rowsID;
    isXlbthisRow = sum(isXlabel(thisRow))>0;
        
    %if Xlabel in this row, Xlabels this row's axes 
    if isXlbthisRow==1
        therows = setdiff(find(thisRow==1),find(isXlabel==1));
        for j = 1 : numel(therows)
            thisAxis = allAxesInFigure(therows(j));
            xlabel(thisAxis,'.','color',...
                get(gcf,'color'))
        end
    end
end

%get coordinates of each axis
[~,~,col] = unique(pos(:,1),'stable');
[~,~,row] = unique(pos(:,2),'stable');

%dead space as percent width and height
deadSpaceRows = (1/nRows)*percDeadSpace;
deadSpaceCols = (1/nCols)*percDeadSpace;

%case dead space covers all the figure width or height
if (nCols+1)*deadSpaceCols >= 1
    maxPossibleDeadSpace = 1/(nCols+1);
    fprintf('%s \n',['With this number of colums and row axes ',...
        'deadSpace should be < ', ...
        num2str(maxPossibleDeadSpace)])
    keyboard
end

%total deadspace
fixSpaceForLabel = 0.035;

%row
totalDeadSpaceRows = (nRows+1)*deadSpaceRows + fixSpaceForLabel;
SpaceforAxesRows = 1 - totalDeadSpaceRows; 

%cols
totalDeadSpaceCols = (nCols+1)*deadSpaceCols + fixSpaceForLabel; 
SpaceforAxesCols = 1 - totalDeadSpaceCols; 

%axes width
widthAxes = SpaceforAxesCols/nCols;
if widthAxes*nCols > 1
    fprintf('%s \n','Some axes are out of the figure')
end

%axes height
heightAxes = SpaceforAxesRows/nRows;
if heightAxes*nRows > 1
    fprintf('%s \n','Some axes are out of the figure')
end

%new axes coordinates
%associate each axis with a position in a template.
coltemp = repmat(1:nCols,1,nRows);
coltemp = coltemp(:);
rowtemp = repmat(1:nRows,nCols,1);
rowtemp = rowtemp(:);

%axes origins
%column
originCols = nan(1,nCols);
for i = 1 : nCols
    originCols(i) = (i-1)*widthAxes + i*deadSpaceCols + fixSpaceForLabel;
end

%rows
originRows = nan(nRows,1);
for i = 1 : nRows
    originRows(i) = (i-1)*heightAxes + i*deadSpaceRows + fixSpaceForLabel;
end
originRows = flipud(originRows);

%axes
XoriginAxes = originCols(ones(nRows,1),:);
XoriginAxes = XoriginAxes';
YoriginAxes = originRows(:,ones(1,nCols));
YoriginAxes = YoriginAxes';

%axes width and heights
widthAxesAll = widthAxes(ones(numAxes,1));
heightAxesAll = heightAxes(ones(numAxes,1));  

%make space for xlabels
heightAxesAll = heightAxesAll - 0.1*heightAxesAll;
YoriginAxes(:) = YoriginAxes(:) + 0.1*YoriginAxes(:);

%axes outerposition
AxesOuterPositions = [XoriginAxes(:) YoriginAxes(:) widthAxesAll heightAxesAll];

%re-position them without dead space in the matched template
%[left bottom width height]
for i = 1 : numAxes
    idx = (SLmakeColumn(col==coltemp(i)) & SLmakeColumn(row==rowtemp(i)));
    set(allAxesInFigure(idx),'position',AxesOuterPositions(i,:));
end







