
% slinitFilename


%slinitFilename('rename')

function  slinitFilename(input)

taskNum = 2;
VarToAddToname = {'myRandomCoh','myMean','myStrength'};

%get stimfiles
stimfiles = dir('*.mat');

%check that stimfiles have not yet been renamed
%they should contain "stim" string
if ~isempty(stimfiles)
    
    if isempty(strfind(stimfiles(1).name,'stim'))
        fprintf('%s \n','stimfiles do not contain "stim". They have probably been renamed already....')
        keyboard
    end
else
    fprintf('%s \n','There are no stimfiles in the directory')
    keyboard
end

%loop over stimfiles and rename them
for i = 1 : length(stimfiles)
    
    %load    
    load(stimfiles(i).name);
    
    %case rename
    if strcmp(input,'rename')
        
        %get parameters
        e = getTaskParameters(myscreen,task);
        
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
            newname = ['s' sprintf('%02d',i) '_' newname VarToAddToname{j} '_' a];
        end
        newname = [newname '.mat'];
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
