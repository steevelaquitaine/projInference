
%slmakeSwitchingVar.m
%
%
%    Author: steeve laquitaine
%      date: 151217
%    status: Complete
%   purpose: classifies estimates in a "switching" variable 
%            Classification boundary is set at equal vector distance between
%            statistical prior mean and displayed motion directions such that 
%           
%            - Estimates near prior/dir are labelled 1/2. Those trials do
%            not contain motion directions at prior mean because direction
%            is prior mean thus classifying estimates as switching-to-prior
%            mean or direction is impossible.
%            - missing trials are set to 0
%            - "mySwitch" when motion direction is at prior (e.g., 225) are
%            labelled -1;
%
%   usage:
%
%     o.sessPath{1}  = '~/data/datafMRI/sltaskdotdirfmri05/s032720160326'; 
%     o.sessPath{2}  = '~/data/datafMRI/sltaskdotdirfmri05/s032720160328'; 
%          o.myGroup = 'Concatenation';
%          o.taskNum = 2;
%        o.priormean = 225;
%
%           o = slmakeSwitchingVar(o)
%
%  Inputs:
%           o.sessPath :  fmri session path: where your fmri data are saved
%                         nifti files etc...
%            o.myGroup :  One of the group for mrtools: 'Raw', 'MotionComp', 
%                         'Concatenation', 'Average'
%            o.taskNum :  1,2 : the task number (see mgl, http://gru.stanford.edu)
%          o.priormean :  135, 225,NaN: experimental prior (NaN for uniform prior)
%
%Note:
%
%   Requires mgl and mrTools distribution to work 
%   (see https://github.com/justingardner)


function [o,behBySess] = slmakeSwitchingVar(o)

mrQuit
fprintf('%s \n','==============================')
fprintf('%s \n',['Prior mean is ' num2str(o.priormean) ' deg'])
fprintf('%s \n','==============================')

tic
for sess = 1 : length(o.sessPath)
    
    %session    
    fprintf('%s%s \n','SESSION:',o.sessPath{sess})
    fprintf('%s \n','==============================')
    o.myPath = o.sessPath{sess};

    %Set concatenated scan
    cd(o.myPath)
    v = mrLoadRet([]);
    v = viewSet(getMLRView,'curGroup',o.myGroup);
    
    %---------------------------------------------------------------
    %loop over stimfiles and get data and task conditions.
    %Update them with new variable (switching and missing responses)
    %---------------------------------------------------------------
    myStimfiles = viewGet(v,'stimfile');
    
    fprintf('%s%s%s \n','(slmakeSwitchingVar) Stimfiles found for "',o.myGroup,'" are:')
    for i = 1: length(myStimfiles)
        [pthi stimFilenm] = fileparts(myStimfiles{i}.filename);
        fprintf('%s \n',stimFilenm)
    end
    
    nStimf = length(myStimfiles);
    fprintf('%s %i %s \n','(slmakeSwitchingVar)',nStimf, 'stim files have been found.')
    fprintf('\n %s \n','--------(slmakeSwitchingVar: add new variables to stimfiles)---------------------')
    
    %init
    o.mySwitch = [];
    o.estimates = [];
    o.coh = [];
    o.directions = [];
    o.dirAtPrior =[];
    
    for stimi = 1 :  nStimf
        
        %stim file
        [fpath,fnm] = fileparts(myStimfiles{stimi}.filename);
        fprintf('%s \n',['(slmakeSwitchingVar) ' fnm])
        fprintf('%s \n','==================================')
        
        %%%--------------------------------------
        %Create "Switching" classification variable
        %%%--------------------------------------
        %load beh. data
        %create "Switching" variable
        % data = loadScan(v);%load behavior
        e = getTaskParameters(myStimfiles{stimi}.myscreen,myStimfiles{stimi}.task);
        e = e{o.taskNum};
        [~,estimatesDeg] = SLcart2polar(cell2mat(e.randVars.prodcoor')); %estimates
        estimatesDeg = round(estimatesDeg);
        
        %create variable indicating when subject responded
        fprintf('%s \n','Setting missing response trials to 0 in variable missing.....')
        missing         = find(isnan(estimatesDeg));
        isresponse      = ones(length(estimatesDeg),1);
        isresponse(missing) = 0;
        
        %label trials as "switched-to-prior" (1) and "switched-to-motion dir" (2)
        motiondir  = e.randVars.myRandomDir;  					          %direction
        %dist est to prior
        d2prior = SLvectors2signedAngle(estimatesDeg,o.priormean,'polar');
        %dist est to dir
        d2dir = SLvectors2signedAngle(estimatesDeg,motiondir,'polar');
        %abs distance of estimate to prior (col 1) and dir (col 2)
        distPandL = abs([d2prior d2dir]);
        mySwitch = zeros(length(estimatesDeg),1);
        %estimates near prior/dir are switch 1/2
        for i = 1 : length(estimatesDeg)
            [~,mySwitch(i,1)] = min(distPandL(i,:));%switch var
        end
        %case uniform prior all switch at motion direction
        if o.priormean == NaN
            mySwitch = ones(1,length(estimatesDeg))*2;
        end
        mySwitch(missing) = 0;%missing
        mySwitch(motiondir==225) = -1; %remove trials with direction at prior        
        nMiss = sum(mySwitch==0);
        nTrials = length(mySwitch);
        fprintf('%s %i %s %i %s %i %s \n','(slmakeSwitchingVar) ',round(nMiss/nTrials*100),'% (' ,nMiss,'/',nTrials,') of the responses are missing.')
        
        %if raw Stim has not been backed up yet create a backup
        cd(o.myPath)
        cd Etc
        %stim file name and path
        [pthi stimFilenm] = fileparts(myStimfiles{stimi}.filename);
        if ~exist([pthi '/rawStimFiles/' stimFilenm '.mat'],'file')
            mkdir rawStimFiles
            copyfile(myStimfiles{stimi}.filename,'rawStimFiles/') %backup raw file
            %cd ..
            fprintf('%s %s %s \n','(slmakeSwitchingVar)',stimFilenm, 'has been backed up in rawStimFiles.')
        end
        
        %------------------------------
        %ADD NEW VARIABLES TO STIM FILE
        %------------------------------
        %add variables to new stimfile (for fMRI analyses) ; create folder and backup raw stimfiles !!!
        %add "isresponse" and "switch" variables
        load(stimFilenm)
        task{o.taskNum}{1}.randVars.isresponse = isresponse;
        task{o.taskNum}{1}.randVars.mySwitch   = mySwitch;
        task{o.taskNum}{1}.randVars.myRandomDirAtPrior = mySwitch==-1;
        task{o.taskNum}{1}.randVars.estimatesDeg = estimatesDeg;
        
        nAddedVar   = 4;
        addedVarnm  = {'isresponse','mySwitch','myRandomDirAtPrior','estimatesDeg'};
        addedVarLen = [nTrials nTrials nTrials nTrials];
        
        for i = 1 : nAddedVar
            task{o.taskNum}{1}.randVars.names_(task{o.taskNum}{1}.randVars.n_ + i) = addedVarnm(i);
            task{o.taskNum}{1}.randVars.varlen_(task{o.taskNum}{1}.randVars.n_ + i) = addedVarLen(i);
        end
        
        task{o.taskNum}{1}.randVars.n_   = task{o.taskNum}{1}.randVars.n_ + nAddedVar;
        save(stimFilenm,'fixStimulus','myscreen','stimulus','task')
        load(stimFilenm)
        fprintf('%s \n','(slmakeSwitchingVar) Adding variables "isresponse", "mySwitch" and "myRandomDirAtPrior" and "estimatesDeg" to the stimfiles')
        
        %output
        o.mySwitch   = [o.mySwitch;task{o.taskNum}{1}.randVars.mySwitch];
        o.estimates  = [o.estimates;estimatesDeg];
        o.dirAtPrior = [o.dirAtPrior;task{o.taskNum}{1}.randVars.myRandomDirAtPrior];
        o.directions = [o.directions;motiondir'];
        o.coh        = [o.coh;e.randVars.myRandomCoh'];  %direction;
        clear e
    end
    
    %sanity check that all variables
    %have same number of rows
    fprintf('%s \n \n','----------- Adding variables done ---------------')
    
    fprintf('%s \n','----------- Checking  ---------------')
    n = unique([length(o.mySwitch) ;length(o.estimates); length(o.dirAtPrior);length(o.directions); length(o.coh)]);
    if length(n)==1
        fprintf('%s \n','(slmakeSwitchingVar) Check: all variables have same # of trials. Good.')
    end
    
    fprintf('%s \n \n','----------- Check done ---------------')
    cd(o.rootpath)
    mrQuit

    %save data each session
    behBySess{sess} = o;
end
toc