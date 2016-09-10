
%slworkflowRenameBehData.m
%
%
%  author: steeve laquitaine
% purpose: rename all data files with their associated condition
%
%   usage:
%
%       slworkflowRenameBehData('s024','sub02',filesBySession,'p225')


function slworkflowRenameBehData(subject,ownSubId,filesBySession,prior)

%% organize data by conditions in separate folders
%move to subject folder
cd(['~/Desktop/sltaskDotDirfMRI05/' subject])
fprintf('%s \n',['(workflow) moved to ~/Desktop/sltaskDotDirfMRI05/' subject])
if isempty(dir(['*' prior]))    
    mkdir(prior)
end
fprintf('%s \n',['(workflow) moved to folders "' prior '"'])

%% rename files by session, run, date and task conditions for analysis
%label files with info for analysis
filesBySessionall = [filesBySession{:}];
%loop over files rename and move them to their prior folder
for file_i = 1 : length(filesBySessionall);
        
    %which file
    filename = filesBySessionall{file_i};
    
    %check file in path and "corrupted" folder
    if isempty(dir(filename))
        if ~isempty(dir(['corrupted/' filename]))            
            fprintf('%s \n',['(slworkflowRenameBehData) ' filename ' is corrupted and was moved to the "corrupted" folder'])
            continue
        else
            fprintf('%s \n',['(slworkflowRenameBehData) ' filename ' is missing'])
            keyboard
        end
    end
    
    %check that file is ok if not go next iteration
    load(filename)
    datai = getTaskParameters(myscreen,task);%speed consuming
    if ~isfield(datai{2},'randVars');
        mkdir corrupted
        movefile(filename,'corrupted')
        fprintf('%s \n',['(slworkflowRenameBehData) skipping' filename ' because corrupted'])
        fprintf('%s \n',['(slworkflowRenameBehData)' filename 'was moved to folder "corrupted"'])
        continue
    end        
    %which session
    for i = 1 : length(filesBySession)
        for j = 1 : length(filesBySession{i})
            if any(strcmp(filename,filesBySession{i}{j})); ses = i; end
        end
    end
    ses = sladdZerobeforeNum2str(ses);
    ses = ['session=' ses];
    %which run
    run = find(strcmp(filename,filesBySessionall));
    run = ['run=' sladdZerobeforeNum2str(run)];
    %which pstd
    load(filename);
    pstd = ['pstd=' num2str(task{2}{1}.randVars.myStrength)];
    %which prior mean
    pmean = ['pmean=' num2str(task{2}{1}.randVars.myMean)];
    %which coherences
    coh = unique(task{2}{1}.randVars.myRandomCoh);
    tm ='coh=[';
    for i = 1 : length(coh)
        tm = [tm num2str(coh(i)) ','];
    end
    tm(end)=']'; tm(tm=='.')=[];
    %which exp
    exp = 'exp=outscanner';
    %which date
    dat = ['date=' filename(1:6)];
    %label file with info
    filepath = slrenameFileForAnalysis(filename,ownSubId,ses,run,pstd,pmean,coh,exp,dat);
    movefile(filepath,prior)
    fprintf('%s \n',['(workflow) moved labelled file to ' prior ' folder'])
end
