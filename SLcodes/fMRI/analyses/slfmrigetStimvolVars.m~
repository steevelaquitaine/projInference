
%slfmrigetStimvolVars.m
%
%
% author: steeve laquitaine
%   date: 160112
% status: complete and tested
%purpose: retrieve variables associated with instances using stimvols 
%         (from getStimvol)
%
%usage:
%
%e.g., 
%     
%     %set params
%     %o.thisSessPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
%     o.thisSessPath = '/Volumes/Transcend/data/sltaskdotdirfmri05/s02520150814/';
%     %o.myGroup  = 'Concatenation';
%     o.myGroup  = 'MotionComp';
%     o.myScan     = 1;
%     o.taskNum    = 2;
%     o.phaseNum   = 1;
%     o.segmentNum = 2;
%     variables = {'myRandomDir','myRandomCoh'};
%       
%     %get stimvols you want variables of
%     cd(o.thisSessPath)
%     v = mrLoadRet([]);
%     v = viewSet(v,'curGroup',o.myGroup,'curScan',o.myScan);;
%     stimvols = getStimvol(v,'myRandomCoh','taskNum',o.taskNum,'phaseNum',...
%           o.phaseNum,'segmentNum',o.segmentNum);
%       
%     %get variable values associated with stimvols
%     [vars,o,db] = slfmrigetStimvolVars(instances,stimvols,variables,o)
% 
%
%inputs: 
%     stimvols: input stimvols
%    variables: variable values we want for stimvols
%            o: fmri session parameters
%
%output:           
%      db: sorted variables and stimvols
%       o: parameters

function [db,o] = slfmrigetStimvolVars(stimvols,variables,o)

%open session
mrQuit
cd(o.thisSessPath)
v = mrLoadRet([]);

%set group & scan
v = viewSet(v,'curGroup',o.myGroup);
v = viewSet(v,'curScan',o.myScan);

%keep
otm = o;

%get all stimvols and variables
fprintf('%s \n','--- Getting all stimvols for variable values found in session scan --')
[d,o.variables] = getStimvol(v,variables,'taskNum',o.taskNum,'phaseNum',...
    o.phaseNum,'segmentNum',o.segmentNum);

%update parameters in o
otm.variables = o.variables;
o = otm;

%make database of all stimvols and associated variables 
%found in the session group and scan for each repeat (trials)
nCond = length(o.variables);
j = zeros(1,nCond);
for ivar = 1 : length(variables)
    
    %init variable values
    dbtmp.(variables{ivar}) = [];
    
    %get variable in variable list
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
        dbtmp.(variables{ivar}) = [dbtmp.(variables{ivar}); repmat(valThisVar,nTrialsThisVar,1)];
    end
    
    %get and sort stimvols matched with associated variables
    dbtmp.stimvols(ivar,:) = cell2mat(d(thisVarpos));
    %[db.stimvols,sorted] = sort(db.stimvols);
    %db.(variables{ivar}) = db.(variables{ivar});    
end

%map stimvols to variables in a database
for i = 1 : length(variables)
    [db.stimvols(i,:),pos(i,:)] = sort(dbtmp.stimvols(i,:));
    db.(variables{i}) = dbtmp.(variables{i})(pos(i,:)); 
end

%check
check = sum(bsxfun(@minus,db.stimvols(end,:),db.stimvols(end,:)));
if check ~= 0
    fprintf('Something wrong with your stimvols mapping to variables')
    keyboard
end
db.stimvols = db.stimvols(1,:);

%get variables for input stimvols
stimvols = sort(cell2mat(stimvols));
[db.stimvols,isort] = intersect(db.stimvols,stimvols);

%map stimvols to variables in a database
for i = 1 : length(variables)
    db.(variables{i}) = db.(variables{i})(isort); 
end

mrQuit






