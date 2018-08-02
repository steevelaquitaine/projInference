
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
%optinos :
%       'sortRunbySub':
%                'run':
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
        
        %case all files are to load
        datadir = dir('*data*.mat');

        %case runs are specified, load them
        if any(strcmp(varargin,'run'))

            %get all files
            a = {datadir.name};
            
            %select input sessions and runs
            sessions = varargin{find(strcmp(varargin,'session'))+1};
            runs = varargin{find(strcmp(varargin,'run'))+1};           
            runsSel = nan(length(runs),1);
            
            %select sessions and runs specified in inputs
            %for this subject
            if any(strcmp(varargin,'sortRunbySub'))
                sort4thisSub = varargin{find(strcmp(varargin,'sortRunbySub'))+1};
                posSub_i = strcmp(subjects{j},sort4thisSub)                    ;
                runs = runs(posSub_i);
                sessions = sessions(posSub_i);
            end
            
            %load these runs for this subject
            for i = 1 : length(runs)                              
                b = find(cellfun(@isempty,strfind(a,sessions{i}))==0);
                c = find(cellfun(@isempty,strfind(a,runs{i}))==0);
                runsSel(i) = intersect(b,c);
                datadirtmp(i) = datadir(runsSel(i));
            end            
            datadir = datadirtmp;
            clear datadirtmp;
        end     
                           
        %check data are there
        if isempty(datadir)
            disp(strcat('No data found in : ',Folds.name(j)))
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
                     
        
        %get each file task variables
        for i = 1 : numel(datalisting)
            
            datalist_i = datalisting(i);
            [est_coor,Trials,filedetails,sub,runid,run_id,runr,priormean,session,Pstd,...
            priormodes,FeatureSample,StimStrength,ReactionTime,estimatedFeaturef4Visu,...
            est_deg,FeatureSample_coord] = slgetTaskVars(datalist_i,...
                est_coor,Trials,filedetails,sub,runid,run_id,runr,priormean,session,Pstd, priormodes,FeatureSample,...
                StimStrength,ReactionTime,varargin{:});           
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
    
end

%case a fraction of the data is wanted
if any(strcmp(varargin,'keepFraction'))
    fractokeep = varargin{find(strcmp(varargin,'keepFraction'))+1};    
    %sanity check    
    if length(fractokeep) < 2 
        fprintf('%s \n','(SLMakeDatabank) "keepFraction" must be a two-entries vector.')
        keyboard
    end
    
    %get fractions
    fractokeepStart = fractokeep(1);    
    fractokeepEnd = fractokeep(2);
    
    %get first and last trials of collected dataset
    N = length(databank.data);
    if fractokeepStart==0
        trial0 = 1;else trial0 = round(fractokeepStart*N);
    end
    trialEnd = round(fractokeepEnd*N);
        
    %collect fraction of the dataset
    databank.data                = databank.data(trial0:trialEnd,:);
    databank.subjects            = databank.subjects(trial0:trialEnd);
    databank.filedetails         = databank.filedetails(trial0:trialEnd);
    databank.session             = databank.session(trial0:trialEnd);
    databank.run                 = databank.run(trial0:trialEnd);
    databank.run_id              = databank.run_id(trial0:trialEnd);
    databank.Trials              = databank.Trials(trial0:trialEnd);
    databank.estimatesDeg        = databank.estimatesDeg(trial0:trialEnd);
    databank.estimatesCoord      = databank.estimatesCoord(trial0:trialEnd);
    databank.Pstd                = databank.Pstd(trial0:trialEnd);
    databank.priormean           = databank.priormean(trial0:trialEnd);
    databank.priormodes          = databank.priormodes(trial0:trialEnd);
    databank.stimFeatureDeg      = databank.stimFeatureDeg(trial0:trialEnd);
    databank.stimFeatureCoord    = databank.stimFeatureCoord(trial0:trialEnd);
    databank.stimStrength        = databank.stimStrength(trial0:trialEnd);
    databank.reactionTime        = databank.reactionTime(trial0:trialEnd);      
    
    fprintf('%s \n','===================================================')
    fprintf('%s \n','(SLMakeDatabank) A fraction of the dataset was set')
    fprintf('%s %i %s %i %s %i %s\n','(SLMakeDatabank) Trials',trial0,'to',trialEnd,'were kept from', N, 'trials')
    
    priors = unique(databank.Pstd);
    npriors = length(priors);
    for i = 1 : npriors; nppriors(i) = sum(databank.Pstd==priors(i));end
    fprintf(['%s ' repmat('%i ',1,npriors) '%s' repmat('%i ',1,npriors) '%s \n'],...
        '(SLMakeDatabank) Priors kepts are:', priors, '(',nppriors,') trials per prior')
    
    cohs = unique(databank.stimStrength);
    ncoh = length(cohs);
    for i = 1 : ncoh; npcoh(i) = sum(databank.stimStrength==cohs(i));end
    fprintf(['%s ' repmat('%.2g ',1,ncoh) '%s' repmat('%i ',1,ncoh) '%s \n'],...
        '(SLMakeDatabank) Cohs kepts are:', cohs, '(',npcoh,') trials per coh')
        
end

%case a stimulus Strength condition to remove is specified, remove all trials of that
%condition
if any(strcmp(varargin,'removeStimStrength'))
    fprintf('%s \n','(SLMakedatabank) removing specified condition')
    strengthToRemove = varargin{find(strcmp(varargin,'removeStimStrength'))+1};
    keepStrengths = (databank.stimStrength ~= strengthToRemove);
else
   keepStrengths = 1:length(databank.stimStrength);
end

%update databank
databank.data                = databank.data(keepStrengths,:);
if isfield(databank,'subjects')
    databank.subjects        = databank.subjects(keepStrengths,:);
end
databank.filedetails         = databank.filedetails(keepStrengths,:);
databank.session             = databank.session(keepStrengths,:);
databank.run                 = databank.run(keepStrengths,:);
if isfield(databank,'run_id')
    databank.run_id          = databank.run_id(keepStrengths,:);
end
databank.Trials              = databank.Trials(keepStrengths,:);
databank.estimatesDeg        = databank.estimatesDeg(keepStrengths,:);
databank.estimatesCoord      = databank.estimatesCoord(keepStrengths,:);
databank.Pstd                = databank.Pstd(keepStrengths,:);
databank.priormean           = databank.priormean(keepStrengths,:);
databank.priormodes          = databank.priormodes(keepStrengths,:);
databank.stimFeatureDeg      = databank.stimFeatureDeg(keepStrengths,:);
databank.stimFeatureCoord    = databank.stimFeatureCoord(keepStrengths,:);
databank.stimStrength        = databank.stimStrength(keepStrengths,:);
databank.reactionTime        = databank.reactionTime(keepStrengths,:);

elapsed = toc(t1);
disp(['elasped:' num2str(elapsed) 'sec'])

%backup the databank
save('datbank','databank');
