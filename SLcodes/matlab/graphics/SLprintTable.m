

%SLprintTable.m
%
% author: Steeve laquitaine
%   date: 140905
%
%purpose: print table
%
%  usage:
%
%   DataMatrix   = [1.5 1.764 2.33; 3.523 0.2 1.34; 1 2 3];
%   rowLabels    = {'row 1', 'row 2', 'row 3'};
%   columnLabels = {'col 1', 'col 2', 'col 3'};
%   SLprintTable(DataMatrix,rowLabels,columnLabels)

function SLprintTable(DataMatrix,rowLabels,columnLabels)

h = figure('units','normalized','Position',[0.3 0.3 .5 .5],...
           'numbertitle','off');       
       
numcol = size(DataMatrix,2);
numrows = size(DataMatrix,2);
scSz = get(0,'ScreenSize');
Width = 0.5*scSz(3);
colwidth =  Width/numcol - 0.1*Width/numcol; 

%table
uitable(h,'Units','normalized','Position',[0 0 1 1],...
              'Data',DataMatrix,... 
              'ColumnName',columnLabels,'RowName',rowLabels,...
              'ColumnWidth',repmat({colwidth},1,numcol),...
              'FontSize',14,...
              'ColumnEditable',repmat(false,1,numcol),...
              'ColumnFormat', {'bank' , 'bank'});
          
          
          