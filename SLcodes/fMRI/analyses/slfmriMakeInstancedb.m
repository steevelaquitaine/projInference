
<<<<<<< HEAD
%usage
=======
%author: steeve laquitaine
%usage:
>>>>>>> 608bdce04de4d0b0a6a30ef97e58ba3016b21d1b
%
%     o.sessPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
%     o.taskNum = 2;
%     o.phaseNum = 1;
%     o.segmentNum = 2;
%     variables = {'myRandomDir','myRandomCoh'};
%     slfmriMakeInstancedb(o,variables)


function slfmriMakeInstancedb(o,variables)

<<<<<<< HEAD
keyboard

=======
>>>>>>> 608bdce04de4d0b0a6a30ef97e58ba3016b21d1b
%open session
cd(o.sessPath)
v = mrLoadRet([]);

%get all stimvols and variables
[d,o.variables] = getStimvol(v,variables,'taskNum',o.taskNum,'phaseNum',...
    o.phaseNum,'segmentNum',o.segmentNum);

<<<<<<< HEAD
dbstack
keyboard
=======
nCond = length(o.variables);
j = zeros(1,nCond);
for ivar = 1 : length(variables)
    
    %init variable values
    db.(variables{ivar}) = [];
    
    %get variable
    a = regexp(o.variables,variables{ivar});
    
    %find variable values    
    thisVarpos = find(~cellfun('isempty',a));    
    
    %get variable values
    nValueThisVar = length(thisVarpos);
    for k = 1 : nValueThisVar
        %get trial #
        valPos = thisVarpos(k);
        volsThisVarValue = d{valPos};
        nTrialsThisVar = length(volsThisVarValue);
        
        %get this variable value
        thisVarValue = o.variables{valPos};
        posValue = regexp(thisVarValue,variables{ivar},'end')+2;
        valThisVar = str2double(thisVarValue(posValue:end));
        
        %add to previous values
        db.(variables{ivar}) = [db.(variables{ivar}); repmat(valThisVar,nTrialsThisVar,1)];
    end  
end



dbstack
keyboard
>>>>>>> 608bdce04de4d0b0a6a30ef97e58ba3016b21d1b
