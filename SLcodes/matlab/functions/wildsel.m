
function idx = wildsel(array)
% WILDSEL - wild card selection from cell array of strings using a GUI.
%   It allows the use of wildcards '*' and '?' and displays only the matching 
%   elements of the cell array.
%   The user selects the entry he/she wants and the index of this
%   element in the original array is returned.
%   The '*' wildcard stands for any number (including zero) of characters.
%   The '?' wildcard stands for a single character.
%
%   Usage:
%       IDX = WILDSEL(ARRAY) returns the index, IDX, of the selected item in the
%       cell array of strings ARRAY.
%
%   Example:
%       A = {'Hello world!'; 'Goodbye world!'; 'Goodbye everyone'};
%       idx = wildsel(A);
%
%   brings up the GUI.  Type the search string in the top edit box, e.g. *world?
%   and press [Return] to display the reduced list.  Click on the item to select
%   and the index of this string in the cell array, A, is returned.

% $Revision: 1.2 $ $Date: 2007/03/01 09:42:47 $ $Author: RIS $
% Richard Stephens (ristephens@theiet.org)

% Make sure the index is always defined - even if the figure is destroyed
idx = 1;

defaultBackground = get(0,'defaultUicontrolBackgroundColor');

hfig = figure('Color',defaultBackground,...
    'Units','characters',...
    'Position',[10,10,50,40],...
    'NumberTitle','off',...
    'Name','Wildcard selection',...
    'Toolbar','none',...
    'MenuBar','none',...
    'WindowStyle','modal',...
    'Tag','wildselFig');

uicontrol('Style','edit',...
    'Units','characters',...
    'Position',[1,35,48,2],...
    'Callback',@search_callback,...
    'Tag','searchStr');
uicontrol('Style','listbox',...
    'Units','characters',...
    'Position',[1,1,48,33],...
    'String',array,...
    'Callback',@list_callback,...
    'Tag','listBox');
uicontrol('Style','text',...
    'Units','characters',...
    'Position',[1,37,48,2],...
    'String','Type search string using wildcards * and ?')

% On creation of the GUI save the entry data
dd.array = array;
dd.iMatch = (1:size(array,1))';
dd.idx = 1;
guidata(hfig, dd);

% Make sure that the current object is the edit box so that
% the user can just start typing
hh = guihandles(hfig);
uicontrol(hh.searchStr);

% Now wait for a selection
uiwait(hfig)

% If the figure has been destroyed then the handle will have disappeared
if ishandle(hfig)
    dd = guidata(hfig);
    idx = dd.idx;
    close(hfig)
end

function search_callback(varargin)
% Re-display the list box strings using only those entries that match
% the selection.
hh = guihandles(gcbo);
dd = guidata(gcbo);

expStr = get(hh.searchStr,'String');

% Replace * and ? with regular expression equivalents
regStr = ['^',strrep(strrep(expStr,'?','.'),'*','.{0,}'),'$'];

% Find the cells that match the regular expression
starts = regexpi(dd.array, regStr);
iMatch = ~cellfun(@isempty, starts);

newList = dd.array(iMatch);
set(hh.listBox, 'String',newList);

% Save results
dd.iMatch = iMatch;
guidata(gcbo,dd);

function list_callback(varargin)
% When the listbox item is selected, save away
% the index from the original array
hh = guihandles(gcbo);
dd = guidata(gcbo);

ilist = get(hh.listBox,'Value');

% We've selected a value from the reduced list, so
% find the corresponding entry in the full list.
idxMatch = find(dd.iMatch);
dd.idx = idxMatch(ilist);

guidata(gcbo,dd);

% and resume the main GUI to return this value
uiresume(hh.wildselFig)