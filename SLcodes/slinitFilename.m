
% slinitFilename
%
% author: steeve laquitaine
%   date: 160420
%purpose: rename stimfiles with task parameters
%         or revert to their original names
% usage: 
%
%   slinitFilename('rename','sub02','sess01',2,{'myRandomCoh','myMean','myStrength'})
%   slinitFilename('revertToOriginal','sub02','sess01',2,{'myRandomCoh','myMean','myStrength'})
%
%
%INPUTS: 
%input : 
%              'rename': rename stimfile with variable parameters
%   ' revertToOriginal': revert stimfile name to original name
%   
%   sub : 'sub01'

function  slinitFilename(input,sub,sess,taskNum,VarToAddToname)

%get stimfiles
stimfiles = dir('*.mat');

%check that stimfiles have not yet been renamed
%they should contain "stim" string
if isempty(stimfiles)        
    fprintf('%s \n','There are no stimfiles in the directory')
    keyboard
end

%create subject folder
%if not done yet
[~,curFolder] = fileparts(pwd);
if ~strcmp(curFolder,sub)
    mkdir(['data/' sub])
end

%loop over stimfiles and rename them
for i = 1 : length(stimfiles)   

    %load    
    load(stimfiles(i).name);
    
    %case rename
    if strcmp(input,'rename')        
        %get parameters
        e = getTaskParameters(myscreen,task);        
        
        %automatically assign a run to the stimfile
        run = ['run' sprintf('%02d',i)];
        
        %write variable info
        newname = [];
        for j = 1 : length(VarToAddToname)
            %get rid of space in variable name
            a = num2str(unique(e{taskNum}.randVars.(VarToAddToname{j})));
            empties = strfind(a,' ');
            if ~isempty(empties)
                a(empties(1))='_';
                a(empties(2:end))=[];
            end
            newname = [newname VarToAddToname{j} '_' a];
        end
        newname = ['data' sprintf('%02d',i) '_' sub '_' sess '_' run '_' newname '.mat'];
        movefile(stimfiles(i).name,newname)
        fprintf('%s \n',['(slinitFilename) rename file "' stimfiles(i).name ...
            '" to "' newname '"'])
        clear newname
    end
    
    %case revert to original name
    %write variable info
    %case rename
    if strcmp(input,'revertToOriginal')
        [~,origname] = fileparts(myscreen.stimfile);
        movefile(stimfiles(i).name,[origname '.mat'])
        fprintf('%s \n',['(slinitFilename) rename file "' stimfiles(i).name ...
            '" to "' origname '"'])
        clear origname
    end
end

%move file to subject folder
%if not done yet
[~,curFolder] = fileparts(pwd);
if ~strcmp(curFolder,sub)
    movefile('*.mat',['data/' sub '/'])
end
