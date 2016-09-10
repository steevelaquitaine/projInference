
%SLexistFolder.m
%
% author: steeve laquitaine
%   date: 15/04/01
%purpose: check if path exists
%
%usage:
%
%       SLexistFolder('Dropbox')


function output = SLexistFolder(folder)

folders = dir;

for i = 1 : length(folders)
    thisFolder = folders(i).name;
    isFolder(i) = strcmp(folder, thisFolder);
    
    if isFolder(i) == 1
        output = 1;
        fprintf(['(SLexistFolder) I found your folder "',folder,'" \n'])
        return
    end
end

