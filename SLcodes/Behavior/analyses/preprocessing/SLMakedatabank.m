
% SLMakedatabank.m
%
%      author: steeve laquitaine
%        date: 140522 (last update 140522)
%     purpose: collect data and experimental conditions for each trial of
%              experiment in a matrix called databank. The matrix is automatically
%              saved in the current directory.
%
%      usage:
%
%ex1:
%               %detailed subject analyses
%               databank = SLMakedatabank({'sub01'},'dataPath',...
%                   '~/Dropbox/myDropbox/Project_withJustin/data/dataPsychophy/proj01_priorStrength/',...
%                   'experiment','vonMisesPrior');
%
%ex2:
%               %quick analysis on loaded data
%               databank = SLMakedatabank;
%
%
%description: For now this code is specifically used with analysis code but
%             I will make it more general. It should be based on the idea
%             that databank is a matrix that saves
%             in column:
%               - data
%               - their associated conditions in other column vectors.
%             in rows:
%               - trials
%
%
%mandatory inputs:
%
%   'dataPath': set data path
%               e.g., 'dataPath','~/'
%   'experiment','vonMisesPrior' (bimodalPrior)
%
%
%selected variables (faster than all variables)
%'variables', then
%       e.g., {data,'d','coh','pstd','priormodes','estimatedFeature'};
%


function databank = SLMakedatabank(subjects,varargin)

%timing
t1 = tic;

%--------------- Quick check data --------------------
%data are already loaded in the workspace
if nargin==0
    slPrintfStr('SLMakedatabank',' Making database from loaded data...')
    
    %get workspace variables
    slGetVarIntoWS

    %get file name     
    [~,datalisting] = fileparts(myscreen.stimfile);
    
    %get data from 'task.mat' in workspace
    datai = getTaskParameters(myscreen,task);%consuming
        
    %get estimated cartesian coordinates(data)
    est_coor = datai{2}.randVars.prodcoor';   
    
    %estimated cartesian coordinates(data)
    [~,est_deg,~] = SLcart2polar(cell2mat(est_coor));
    
    %convert coordinates to angles (degree)
    coortmp = cell2mat(est_coor);
    estimatedFeatureTmp = SLgetangle(coortmp(:,1),coortmp(:,2));
    estimatedFeatureTmp(estimatedFeatureTmp==0) = 360;
    estimatedFeaturef4Visu = num2cell(estimatedFeatureTmp);
    
    %calculate the number of trials
    numTrials = numel(est_coor);
    Trials = num2cell(1:numTrials)';
    
    %get the file name
    filedetails = repmat({datalisting},numTrials,1);        

    %get the std of the prior from the file nm
    if isempty(isfield(datai{2}.randVars,'myStrength'))
        fprintf('%s /n','(SLMakedatabank) The filename does not contains "myStrength" information...')
        return
    end        
    %check that mystrength is in deg and not in concentration
    if unique(datai{2}.randVars.myStrength) == 80
        Pstd_thisT = 80;
    elseif unique(datai{2}.randVars.myStrength) == 40
        Pstd_thisT = 40;
    elseif unique(datai{2}.randVars.myStrength) == 40
        Pstd_thisT = 20;
    elseif unique(datai{2}.randVars.myStrength) == 40
        Pstd_thisT = 10;
    end
    %case concentration parameter is input    
    if floor(datai{2}.randVars.myStrength*10)/10==0.7
        Pstd_thisT = 80;
    elseif floor(datai{2}.randVars.myStrength*10)/10==2.7
        Pstd_thisT = 40;
    elseif floor(datai{2}.randVars.myStrength*10)/10==8.7
         Pstd_thisT = 20;
    elseif floor(datai{2}.randVars.myStrength*10)/10==33.3
         Pstd_thisT = 10;
    end
    Pstd = num2cell(repmat(Pstd_thisT,numTrials,1));    
    
    %get the prior mean
    if isempty(isfield(datai{2}.randVars,'myMean'))
        fprintf('%s /n','(Makedatabank) the filename does not contains "Prior mean" information...')
        return
    end   
    priormean = num2cell(repmat(unique(datai{2}.randVars.myMean),numTrials,1));
    
    %get the modes of the distribution   
    priormodes = repmat({225},numTrials,1);
    
    %get the Displayed directions (group 1)
    %when design is balanced across stimStrengths (same trial number).
    if isfield(datai{1,2}.parameter,'dir')
        
        FeatureSample = num2cell(datai{1,2}.parameter.dir)';
        
        %when design is adjusted according to stimStrengths (different trial nb).
    elseif isfield(datai{1,2}.randVars,'myRandomDir')
        
        FeatureSample = num2cell(datai{1,2}.randVars.myRandomDir)';       
        
    elseif isfield(datai{1,2}.randVars,'myRandomloc')
        
        FeatureSample = num2cell(datai{1,2}.randVars.myRandomloc)';
        
    elseif isfield(datai{1,2}.randVars,'myRandomFeature')
        
        FeatureSample = num2cell(datai{1,2}.randVars.myRandomFeature)';       
        
    end
    
    %Convert the displayed dirs from degrees to cartesian
    %coordinates
    r = 2.5; %radius of the circular patch
    FeatureSample_coord = nan(size(FeatureSample,1),1);
    FeatureSample_coord = num2cell(FeatureSample_coord);
    for jk=1: size(FeatureSample,1)
        FeatureSample_coord{jk} = polar2cartesian(FeatureSample{jk},r);
    end
    
    %get the sample stimulus strength (group 2)
    %when design is balanced across coherences (same trial number).
    if isfield(datai{1,2}.parameter,'coherence')
        
        StimStrength = num2cell(datai{1,2}.parameter.coherence)';
        
        %when design is adjusted according to stimulus strength (different trial nb).
    elseif isfield(datai{1,2}.randVars,'myRandomCoh')
        
        StimStrength = num2cell(datai{1,2}.randVars.myRandomCoh)';       
        
    elseif isfield(datai{1,2}.randVars,'myRandomCon')
        
        StimStrength = num2cell(datai{1,2}.randVars.myRandomCon)';
        
    elseif isfield(datai{1,2}.randVars,'myRandomStimStrength')
        
        StimStrength = num2cell(datai{1,2}.randVars.myRandomStimStrength)';
        
    end
    
    %Get Reaction times
    if isfield(datai{1,2},'reactionTime')
        
        ReactionTime = num2cell(datai{1,2}.reactionTime)';
        
        %when reaction times were not collected set vector of "nan".
    elseif ~isfield(datai{1,2},'reactionTime')
        
        ReactionTime = nan(numel(FeatureSample),1);
        
    end
   
    
    %create our databank
    %enter the names of the variable. Prgressively replace by structure.
    %Code easier to read.
    %OLD METHOD
    databank.nm=[
        {'filedetails'  },...
        {'run'} ...
        {'session'} ...
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
    runr = repmat({1},numTrials,1);
    session = repmat({1},numTrials,1);

    databank.data = [
        filedetails  ...
        runr ...
        session ...
        Trials(:)       ...
        est_coor(:)     ...
        Pstd(:)         ...
        priormean(:)    ...
        FeatureSample(:)  ...
        StimStrength(:)     ...
        estimatedFeaturef4Visu   ...
        FeatureSample_coord...
        ReactionTime(:) ...
        priormodes];
    
    %NEW METHOD
    %store data over subjects
    databank.data                = databank.data;
    databank.filedetails         = filedetails  ;
    databank.run                 = runr ;
    databank.session             = session  ;
    databank.Trials              = cell2mat(Trials);
    databank.estimatesDeg        = est_deg;
    databank.estimatesCoord      = est_coor;
    databank.Pstd                = cell2mat(Pstd);
    databank.priormean           = cell2mat(priormean);
    databank.priormodes          = priormodes;
    databank.stimFeatureDeg      = cell2mat(FeatureSample);
    databank.stimFeatureCoord    = FeatureSample_coord;
    databank.stimStrength        = cell2mat(StimStrength);
    databank.reactionTime        = cell2mat(ReactionTime);
    
    %get rid of missing data
    missingData = isnan(databank.estimatesDeg);
    
    %update databank
    databank.data                = databank.data(~missingData,:);
    databank.filedetails         = databank.filedetails(~missingData,:);
    databank.run                 = databank.run(~missingData,:);
    databank.session             = databank.session(~missingData,:);
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
    
    %backup the databank
    save('datbank','databank');        
    slPrintfStr('SLMakedatabank','..done.')
else

    
    %------------- subject detailed analyses --------------------
    
    
    %status
    fprintf('\n (SLMakedatabank) I am creating a database... \n')
    
    %get what's inside varargin
    if length(varargin)==1
        varargin = varargin{1};
    end        
    
    %get data path
    if any(strcmp(varargin,'dataPath'))
        dataPath = varargin{find(strcmp(varargin,'dataPath'))+1};
        cd(dataPath)
    elseif any(strcmp(varargin,'inpath'))
        %create structured path (/data/sub..)
        mkdir(['data/' subjects{:}])
        %check if mat and edf
        if ~isempty(dir('*.mat')) == 1
            copyfile([pwd '/*.mat'],['data/' subjects{:} '/'])
        end        
        if ~isempty(dir('*.edf')) == 1
            copyfile([pwd '/*.edf'],['data/' subjects{:} '/'])
        end       
        fprintf('%s \n',['(SLMakedatabank) Your files were moved to data/' subjects{:} '/ path for analyses'])
        cd('data/')
        dataPath = pwd;               
        fprintf('%s \n','(SLMakedatabank) You are in data path...')
    else
        fprintf('%s \n','(SLMakedatabank) You need to set dataPath')
        return
    end           
    
    %set the directory to analyse
    %case we are already in data dir, go to parent dir
    a = pwd;
    if strcmp(a(end-4:end),'/data')
        pathFold = pwd;
        
    else
        %otherwise go to /data directory
        pathFold = strcat(pwd,'/data');
    end
    d = dir(pathFold);
    
    %check that directory is ok
    if isempty(d)
        fprintf('(SLMakedatabank) Something wrong with you path. It should look like Exp0.../data/sub0..')
        keyboard
    end
    
    %collect subjects
    datatmp.subjects = subjects;
    
    %existing subjects
    existSub = dir('sub*');
    
    %find the subfolders to analyse
    for ii = 1 : numel(datatmp.subjects)
        
        if isempty(find(strcmp({d.name},datatmp.subjects(ii))))
            warning(['(SLMakedatabank) !!WARNING!! There are no data for ',datatmp.subjects{ii}])
            fprintf(['(SLMakedatabank) Available data are for : ', '\n'])
            
            %list existing subjects
            for i = 1 : length(existSub)
                fprintf([existSub(i).name,' \n'])
            end
            keyboard
        end
        isubType(ii) = find(strcmp({d.name},datatmp.subjects(ii)));
    end
    
    %get the names of the subfolders
    Folds.name = {d(isubType).name}';
    
    %check if correct
    disp(Folds.name)
    %uiwait(msgbox('Check out the command window for the folder names.'));
    
    %count the number of subfolders
    Folds.nb = numel(Folds.name);
    
    %initialize database variables
    databank.data                = [];
    databank.subjects            = [];
    databank.filedetails         = [];
    databank.session             = [];
    databank.run                 = [];
    databank.run_id              = [];    
    databank.Trials              = [];
    databank.estimatesDeg        = [];
    databank.estimatesCoord      = [];
    databank.Pstd                = [];
    databank.priormean           = [];
    databank.priormodes          = [];
    databank.stimFeatureDeg      = [];
    databank.stimFeatureCoord    = [];
    databank.stimStrength        = [];
    databank.reactionTime        = [];    
    
    runid = 0;
    
    %loop over the subfolders to analyse
    for j=1:Folds.nb
        
        %directory of the subfolder to open
        Folds.path(j) = strcat(pathFold,'/',Folds.name(j));
        
        %switch to the subfolder to open
        cd(Folds.path{j})
        
        %delete svn files (if exist)
        !rm ._*;
        
        %specify the files to load
        datadir = dir('*data*.mat');%directory
                
        %check if there are data in the directory
        if isempty(datadir)
            disp(strcat('No data were found in directory: ',Folds.name(j)))
            return
        end
        
        %collect the name of each file
        datatmp.nm = {datadir.name};        %filenms
        
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
        run_id=[];
        
        
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
        end
        
        
        %create our databank
        %enter the names of the variable. Prgressively replace by structure.
        %Code easier to read.
        %OLD METHOD
        databank.nm=[
            {'filedetails'  },...
            {'run'          },...
            {'run_id'},...
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
            num2cell(run_id(:)) ...
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
        databank.subjects            = [databank.subjects       ; sub];
        databank.filedetails         = [databank.filedetails    ; filedetails];
        databank.session             = [databank.session        ; session];
        databank.run                 = [databank.run            ; runr];
        databank.run_id              = [databank.run_id         ; run_id];
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
    end
            
    %get rid of missing data
    missingData = isnan(databank.estimatesDeg);
    
    %update databank
    databank.data                = databank.data(~missingData,:);
    databank.subjects            = databank.subjects(~missingData,:);
    databank.filedetails         = databank.filedetails(~missingData,:);
    databank.session             = databank.session(~missingData,:);
    databank.run                 = databank.run(~missingData,:);
    databank.run_id              = databank.run_id(~missingData,:);
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
    
end
elapsed = toc(t1);
disp(['elasped:' num2str(elapsed) 'sec'])
