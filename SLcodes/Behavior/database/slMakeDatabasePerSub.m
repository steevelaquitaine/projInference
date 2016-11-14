
%slMakeDatabasePerSub.m
%
%
%Working on it .....
%
%

function slMakeDatabasePerSub

%directory of the subfolder to open
Folds.path = strcat(pathFold,'/',Folds.name);

%switch to the subfolder to open
cd(Folds.path)

%delete svn files (if exist)
!rm ._*;

%specify the files to load
datadir = dir('*data*.mat');%directory

%check if there are data in the directory
if isempty(datadir)
    disp(strcat('No data were found in directory: ',Folds.name)
    return
end

%collect the name of each file
datatmp.nm={datadir.name};        %filenms

%count the number of file
datatmp.nb=numel(datatmp.nm);     %num`ber

%set the variables in the databank (a column each)
filedetails={};
sub=[];
Trials={};
session={};
runr={};
est_coor={};
Pstd={};
priormean={};
priormodes={};
priormean={};
FeatureSample={};
FeatureSample_coord={};
StimStrength={};
ReactionTime={};


%check if data have been specified in the function argin
if ~isempty(datatmp.nm)
datalisting = datatmp.nm;
%   datalisting={datalisting};
%   display(["--- A databank is being created with the data specified ----']);

%if data have not been specified,gather data from directory
else
%remove possible svn files
datalisting = dir('*data*.mat');
datalisting = {datalisting.name};
%   display(['--- No data have been specified. A databank is being created with all data in the directory and analysed ---']);
end

%tic

%get each file data
for i = 1 : numel(datalisting)

%load the files and get data from 'task.mat' in the workspace
load(datalisting{i});
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
filedetails_i=repmat(datalisting(i),numTrials,1);
filedetails = [filedetails;filedetails_i];

%get the subject
if isempty(strfind(filedetails_i{1},'sub'))
fprintf('%s /n','(SLMakedatabank) the filename does not contains "subject" information')
return
end
sub_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},...
'sub')+3:strfind(filedetails_i{1},'sub')+4));
sub = [sub;repmat(sub_thisT,numTrials,1)];

%get the run
if isempty(strfind(filedetails_i{1},'run'))
fprintf('%s /n','(SLMakedatabank) the filename does not contains "run" information')
return
end
run_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},...
'run')+3:strfind(filedetails_i{1},'run')+4));
runr = [runr ; num2cell(repmat(run_thisT,numTrials,1))];

%get session
if isempty(strfind(filedetails_i{1},'sess'))
fprintf('%s /n','(SLMakedatabank) the filename does not contains "session" information...')
return
end
session_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},...
'sess')+5));
session = [session ; num2cell(repmat(session_thisT,numTrials,1))];

%get the std of the prior from the file nm;
if isempty(strfind(filedetails_i{1},'Pstd'))
fprintf('%s /n','(SLMakedatabank) the filename does not contains "Pstd" information...')
return
end
Pstd_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2));%std of the prior

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
if isempty(strfind(filedetails_i{1},'mean'))
fprintf('%s /n','(Makedatabank) the filename does not contains "Prior mean" information...')
return
end
priormean_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'mean')+4:strfind(filedetails_i{1},'mean')+6));%run
priormean=[priormean;num2cell(repmat(priormean_thisT,numTrials,1))];

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
end


%create our databank
%enter the names of the variable. Prgressively replace by structure.
%Code easier to read.    
%OLD METHOD
databank.nm=[
{'filedetails'  },...
{'run'          },...
{'session'      },...
{'Trials'       },...
{'est_coor'},...
{'Pstd'         },...
('priormean'    ),...
{'FeatureSample'   },...
{'StimStrength'          },...
{'estimatedFeature'      },...
{'FeatureSample_coor'},...
{'ReactionTime'},...
{'priormodes'}];

%store the data in each column
databanki.data=[
filedetails(:)  ...
runr(:)         ...
session(:)      ...
Trials(:)       ...
est_coor(:)     ...
Pstd(:)         ...
priormean(:)    ...
FeatureSample(:)   ...
StimStrength(:)          ...
estimatedFeaturef4Visu   ...
FeatureSample_coord...
ReactionTime(:) ...
priormodes];


%NEW METHOD
%store data over subjects
databank.data                = [databank.data           ; databanki.data];
databank.subjects            = [databank.subjects       ; sub           ];
databank.filedetails         = [databank.filedetails    ; filedetails   ];
databank.session             = [databank.session        ; session   ];
databank.run                 = [databank.run            ; runr   ];
databank.Trials              = [databank.Trials         ; cell2mat(Trials) ];
databank.estimatesDeg        = [databank.estimatesDeg   ; est_deg  ];
databank.estimatesCoord      = [databank.estimatesCoord ; est_coor   ];
databank.Pstd                = [databank.Pstd           ; cell2mat(Pstd)  ];
databank.priormean           = [databank.priormean      ; cell2mat(priormean) ];
databank.priormodes          = [databank.priormodes     ; priormodes   ];
databank.stimFeatureDeg      = [databank.stimFeatureDeg ; cell2mat(FeatureSample) ];
databank.stimFeatureCoord    = [databank.stimFeatureCoord ; FeatureSample_coord   ];
databank.stimStrength        = [databank.stimStrength   ; cell2mat(StimStrength)   ];
databank.reactionTime        = [databank.reactionTime   ; cell2mat(ReactionTime) ];


%get rid of missing data
missingData = isnan(databank.estimatesDeg);

%update databank
databank.data                = databank.data(~missingData,:);     
databank.subjects            = databank.subjects(~missingData,:);    
databank.filedetails         = databank.filedetails(~missingData,:);    
databank.session             = databank.session(~missingData,:);      
databank.run                 = databank.run(~missingData,:);         
databank.Trials              = databank.Trials(~missingData,:);      
databank.estimatesDeg        = databank.estimatesDeg(~missingData,:);    
databank.estimatesCoord      = databank.estimatesCoord(~missingData,:);    
databank.Pstd                = databank.Pstd(~missingData,:);          
databank.priormean           = databank.priormean(~missingData,:);        
databank.priormodes          = databank.priormodes(~missingData,:);       
databank.stimFeatureDeg      = databank.stimFeatureDeg(~missingData,:);    
databank.stimFeatureCoord    = databank.stimFeatureCoord(~missingData,:);    
databank.stimStrength        = databank.stimStrength(~missingData,:);       
databank.reactionTime        = databank.reactionTime(~missingData,:);    

%make sure data are integers
databank.estimatesDeg = round(databank.estimatesDeg);

%make sure 360 and 0 degrees are same
databank.estimatesDeg(databank.estimatesDeg==0) = 360;


%data are organized and saved in a file called datbank in the directory
%switch back to the mother directory
cd(pathFold)

%backup the databank
save('datbank','databank');