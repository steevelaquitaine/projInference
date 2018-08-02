

function [est_coor,Trials,filedetails,sub,runid,run_id,runr,priormean,session,Pstd,...
    priormodes,FeatureSample,StimStrength,ReactionTime,estimatedFeaturef4Visu,...
    est_deg,FeatureSample_coord] = slgetTaskVars(datalisting,...
    est_coor,Trials,filedetails,sub,runid,run_id,runr,priormean,session,Pstd, priormodes,FeatureSample,...
    StimStrength,ReactionTime,varargin)


%load the files and get data from 'task.mat' in the workspace
% load(datalisting{i});
load(datalisting{1});
datai = getTaskParameters(myscreen,task);%speed consuming

%get the estimated cartesian coordinates(data)
est_coor_i = datai{2}.randVars.prodcoor';
est_coor = [est_coor ; est_coor_i];

%estimated cartesian coordinates(data)
[~,est_deg,~] = SLcart2polar(cell2mat(est_coor));

%convert coordinates to angles (degree)
coortmp=cell2mat(est_coor);
estimatedFeatureTmp = SLgetangle(coortmp(:,1),coortmp(:,2));
estimatedFeatureTmp(estimatedFeatureTmp==0) = 360;
estimatedFeaturef4Visu=num2cell(estimatedFeatureTmp);

%calculate the number of trials
numTrials=numel(est_coor_i);
Trials_i=num2cell(1:numTrials)';
Trials=[Trials;Trials_i];

%get the file name
%filedetails_i=repmat(datalisting(i),numTrials,1);
filedetails_i=repmat(datalisting,numTrials,1);
filedetails = [filedetails;filedetails_i];

%get the subject
if isempty(strfind(filedetails_i{1},'sub'))
    fprintf('%s /n','(SLMakedatabank) the filename does not contains "subject" information')
    return
end
sub_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},...
    'sub')+3:strfind(filedetails_i{1},'sub')+4));
sub = [sub;repmat(sub_thisT,numTrials,1)];

%assign a run id
runid = runid + 1;
run_id = [run_id ; repmat(runid,numTrials,1)];

%get the run
if isempty(strfind(filedetails_i{1},'run'))
    fprintf('%s /n','(SLMakedatabank) the filename does not contains "run" information')
    keyboard
end
run_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},...
    'run')+3:strfind(filedetails_i{1},'run')+4));
runr = [runr ; num2cell(repmat(run_thisT,numTrials,1))];

%get session
if isempty(strfind(filedetails_i{1},'sess'))
    fprintf('%s /n','(SLMakedatabank) the filename does not contains "session" information...')
    keyboard
end
session_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},...
    'sess')+5));
session = [session ; num2cell(repmat(session_thisT,numTrials,1))];



%get the std of the prior from the file nm;
if isempty(strfind(filedetails_i{1},'Pstd')) && isempty(strfind(filedetails_i{1},'myStrength'))
    fprintf('%s /n','(SLMakedatabank) the filename does not contains "Pstd" nor "myStrength" information...')
    keyboard
end
if ~isempty(strfind(filedetails_i{1},'Pstd'))
    Pstd_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2));%std of the prior
elseif ~isempty(strfind(filedetails_i{1},'myStrength'))
    Pstd_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},'myStrength')+1+numel('myStrength'):strfind(filedetails_i{1},'myStrength')+1+numel('myStrength')+2));%std of the prior
end

%get the std of the prior from the file nm when prior is uniform;
if isempty(strfind(filedetails_i{1},'inf'))==0
    Pstd_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2));%std of the prior
end

%check if the Prior's std was found in the name
if isnan(Pstd_thisT)
    fprintf('%s /n','(SLMakedatabank) WARNING ! The code does not find "Pstd" value in the filename...')
    keyboard
end
Pstd = [Pstd;num2cell(repmat(Pstd_thisT,numTrials,1))];

%get the prior's mean
if isempty(strfind(filedetails_i{1},'mean')) && isempty(strfind(filedetails_i{1},'myMean'))
    fprintf('%s /n','(Makedatabank) the filename does not contains "Prior mean" information...')
    keyboard
end
if ~isempty(strfind(filedetails_i{1},'mean'))
    priormean_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},'mean')+4:strfind(filedetails_i{1},'mean')+6));%run
    priormean = [priormean;num2cell(repmat(priormean_thisT,numTrials,1))];
elseif ~isempty(strfind(filedetails_i{1},'myMean'))
    priormean_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},'myMean')+7:strfind(filedetails_i{1},'myMean')+9));%run
    priormean = [priormean;num2cell(repmat(priormean_thisT,numTrials,1))];
end

%get the modes of the distribution (uni/bimodal ?)
%check that experiment type is input
if sum(strcmp(varargin,'experiment'))
    
    %(case von Mises prior experiment)
    %-------------------------------
    if strcmp(varargin(find(strcmp(varargin,'experiment'))+1),'vonMisesPrior')
        priormodes_thisT=225;
        priormodes=[priormodes;repmat({priormodes_thisT'},numTrials,1)];
    end
    
    %(case bimodal prior experiment)
    %-------------------------------
    if strcmp(varargin(find(strcmp(varargin,'experiment'))+1),'bimodalPrior')
        
        %if saved data contain information about prior modes
        if isfield(task{2}{1}.randVars,'myModes')
            %get the modes of the priors used in the experiment
            priormodes_thisT = task{2}{1}.randVars.myModes';
            priormodes = [priormodes;repmat({priormodes_thisT'},numTrials,1)];
        else
            fprintf('% s \n',['(SLMakedatabank) WARNING !! task{2}{1}.randVars',...
                ' does not contain the field "myModes"',...
                '. Are you sure those data correspond to the'...
                ' bimodal experiment ?'])
            keyboard
        end
    end
    
    %if the experiment is the Gaussian prior experiment, default
    if strcmp(varargin(find(strcmp(varargin,'experiment'))+1),'vonMisesPrior')
    end
    
    %if the type of experiment was not input, warn and stop
else
    fprintf('%s \n','(SLMakedatabank) Please indicate the "experiment"...(e.g. "experiment","vonMisesPrior"')
    fprintf('%s \n','(SLMakedatabank) Your choices are: ')
    fprintf('%s \n','(SLMakedatabank)   - vonMisesPrior ')
    fprintf('%s \n','(SLMakedatabank)   - bimodalPrior ')
    keyboard
end

%get the Displayed directions (group 1)
%when design is balanced across stimStrengths (same trial number).
if isfield(datai{1,2}.parameter,'dir')
    
    FeatureSample_i=num2cell(datai{1,2}.parameter.dir)';
    FeatureSample = [FeatureSample;FeatureSample_i];
    
    %when design is adjusted according to stimStrengths (different trial nb).
elseif isfield(datai{1,2}.randVars,'myRandomDir')
    
    FeatureSample_i = num2cell(datai{1,2}.randVars.myRandomDir)';
    FeatureSample = [FeatureSample;FeatureSample_i];
    
elseif isfield(datai{1,2}.randVars,'myRandomloc')
    
    FeatureSample_i = num2cell(datai{1,2}.randVars.myRandomloc)';
    FeatureSample = [FeatureSample;FeatureSample_i];
    
elseif isfield(datai{1,2}.randVars,'myRandomFeature')
    
    FeatureSample_i = num2cell(datai{1,2}.randVars.myRandomFeature)';
    FeatureSample = [FeatureSample;FeatureSample_i];
    
end

%Convert the displayed dirs from degrees to cartesian
%coordinates
r = 2.5; %radius of the circular patch
FeatureSample_coord = nan(size(FeatureSample,1),1);
FeatureSample_coord = num2cell(FeatureSample_coord);
for jk=1: size(FeatureSample,1)
    FeatureSample_coord{jk}=polar2cartesian(FeatureSample{jk},r);
end
%         for jk=1: size(FeatureSample,1)
%             FeatureSample_coord{jk,1}=num2cell(polar2cartesian(FeatureSample{jk},r));
%         end
%         for jk=1:size(FeatureSample,1)
%             FeatureSample_coord(jk)=num2cell(polar2cartesian(FeatureSample{jk},r));
%         end

%get the sample stimulus strength (group 2)
%when design is balanced across coherences (same trial number).
if isfield(datai{1,2}.parameter,'coherence')
    
    StimStrength_i = num2cell(datai{1,2}.parameter.coherence)';
    StimStrength = [StimStrength; StimStrength_i];
    
    %when design is adjusted according to stimulus strength (different trial nb).
elseif isfield(datai{1,2}.randVars,'myRandomCoh')
    
    StimStrength_i = num2cell(datai{1,2}.randVars.myRandomCoh)';
    StimStrength = [StimStrength; StimStrength_i];
    
elseif isfield(datai{1,2}.randVars,'myRandomCon')
    
    StimStrength_i = num2cell(datai{1,2}.randVars.myRandomCon)';
    StimStrength = [StimStrength; StimStrength_i];
    
elseif isfield(datai{1,2}.randVars,'myRandomStimStrength')
    
    StimStrength_i = num2cell(datai{1,2}.randVars.myRandomStimStrength)';
    StimStrength = [StimStrength; StimStrength_i];
    
end

%Get Reaction times
if isfield(datai{1,2},'reactionTime')
    
    ReactionTime_i=num2cell(datai{1,2}.reactionTime)';
    
    %when reaction times were not collected set vector of "nan".
elseif ~isfield(datai{1,2},'reactionTime')
    
    ReactionTime_i=nan(numel(FeatureSample_i),1);
    
end
ReactionTime =[ReactionTime; ReactionTime_i];
