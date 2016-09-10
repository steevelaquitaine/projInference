


%--------------------------------------------------------------------------
% abbreviations used in the code
%--------------------------------------------------------------------------
% see http://www.all-acronyms.com/reverse/LVLS
% FA    : factors
% LVLS  : levels
% T     : trial
% est   : estimate





%--------------------------------------------------------------------------
% Features of the code
%--------------------------------------------------------------------------

% "FAoflvl2del", "lvl2del"
% You can enter the factor(FAoflvl2del) and the corresponding level (lvl2del)
% that you want to discard from the analysis.


% Data are linearized by adding or substracting 360 if needed. Thus
% modeling is possible.





%--------------------------------------------------------------------------
% To do
%--------------------------------------------------------------------------
% clean the code





%--------------------------------------------------------------------------
% The models used in the code
%--------------------------------------------------------------------------

% Currently we made the assumption that the noise in the sensory evidence (std(LLH))
% is represented linearly.
% (see  Hurlimann, F., Kiper, D. C., & Carandini, M. (2002). 
% Testing the Bayesian model of perceived speed. Vision Res, 42(19), 2253?2257.). 





% Different assumptions are made concerning the representation of the noise
% in the sensory evidence (std(LLH))

% #Assumption 1 - model 1.
% width LLH = Sreal * Co^n/(Co^n + C50^n) with n=2 if the representation of
% coherence is assumed to be linear

% #Assumption 2 - model 2 (simplest)
% width LLH = Sreal * C

% #Assumption 2 - model 2 (simplest)
% width LLH = Sreal * C + So (So ~= 0)







function [ minSSE,modelPred,freeParamsBkp, freePara2, Rsquared, SSE_bkp, k_bkp] = modelBayesInf2(data, FAs, fig,FAoflvl2del,lvl2del )
%BAYESINFERENCE Summary of this function goes here
%   Detailed explanation goes here

% Input:
% - coherence (c)
% - prior's mean (uPr) and standard deviation (sPr)
% - displayed directions

% Output:
% - estimated directions.

% Create a databank
Makedatabank(data,FAoflvl2del,lvl2del);
[inputLSfitting, FA] = initFAs(FAs,fig);

% Set the initial free parameters for the LLH
alpha0 = 0;%: 50: 1000;
beta0  = 1;%: 50: 100;
n0     = 2;%-100: 0: 100;
freePara2 = [alpha0 beta0 n0];
% note" I tried many initial parameters. These are the only one that allow
% convergence in a reasonable amount of time.

% Set the free parameter for the bayesian operation
% freeParams =[];
% minSSELLH = [];

% options = optimset('MaxFunEvals',1000,'MaxIter',1000);
freePara2 = fmincon( @(freePara2) LSfitting(inputLSfitting, freePara2 ),...
    [alpha0 beta0 n0],... % initial parameters
    [],[],...             
    [],[],...
    [0    0    -inf],... %lower bound
    [+inf +inf +inf],... % upper bound
    [])%,...
%     options);

%Generate best model's predictions
[minSSELLH, freePara2,minSSE,modelPred, udata, taskParams, freeParamsBkp, Rsquared, SSE_bkp, k_bkp] = LSfitting(inputLSfitting,freePara2);


% [minSSE, modelPred, udata, taskParams, freeParamsBkp, Rsquared, SSE_bkp, k_bkp] = LSfitting(inputLSfitting );
[fig1] = grapDataAndModel(modelPred, udata, taskParams, FA,fig);


%% Create a matrix of cells in which columns are the variables (databank)
function Makedatabank(data,FAoflvl2del,lvl2del)
global databank

% set the directory to analyse
pathFold = pwd;
d = dir(pathFold);

% collect subjects
datatmp.subjects = data.subjects;

% find the subfolders to analyse
for ii = 1 : numel(datatmp.subjects)
    isubType(ii) = find(strcmp({d.name},datatmp.subjects(ii)));
end

% get the names of the subfolders
Folds.name = {d(isubType).name}';

% check if correct
disp(Folds.name)
% uiwait(msgbox('Check out the command window for the folder names.'));

% count the number of subfolders
Folds.nb = numel(Folds.name);

databank.data = [];

%% Loop over the subfolders to analyse
for j = 1: Folds.nb
    
    % directory of the subfolder to open
    Folds.path(j)= strcat(pathFold,'/',Folds.name(j));
    
    % switch to the subfolder to open
    cd(Folds.path{j})

    % delete svn files (if exist)
    !rm ._*;

    % specify the files to load
    datadir = dir('*data*.mat'); % directory
    
    % check if there are data in the directory
    if isempty(datadir)
        disp(strcat('No data were found in directory: ',Folds.name(j)))
        return
    end
    
    % collect the name of each file
    datatmp.nm = {datadir.name};         % filenms
    
    % count the number of file
    datatmp.nb = numel(datatmp.nm);      % num`ber
    
    % set the variables in the databank (a column each)
    filedetails = {};
    Trials = {};
    session = {};
    runr = {};
    est_coor = {};
    Pstd = {};
    priormean = {};
    sample_dir = {};
    coh = {};
    

    % check if data have been specified in the function argin
    if ~isempty(datatmp.nm)
        datalisting = datatmp.nm;
        %     datalisting = {datalisting};
        %     display(["--- A databank is being created with the data specified ----']);
        
        % if data have not been specified, gather data from directory
    else
        % remove possible svn files
        datalisting = dir('*data*.mat');
        datalisting = {datalisting.name};
        %     display(['--- No data have been specified. A databank is being created with all data in the directory and analysed ---']);
    end
    
    % tic
    
    % loop over the files and collect their data
    for i = 1 : numel(datalisting)
        
        % load the files and get data from 'task.mat' in the workspace
        load(datalisting{i});
        datai = getTaskParameters(myscreen,task); %speed consuming
        
        % get the estimated cartesian coordinates (data)
        est_coor_i = datai{2}.randVars.prodcoor';
        est_coor = [est_coor; est_coor_i];
        
        % convert coordinates to angles (degree)
        coortmp = cell2mat(est_coor);
        est_dirf4Visu = num2cell(getangle(coortmp(:,1),coortmp(:,2)));
        
        % calculate the number of trials
        numTrials = numel(est_coor_i);
        Trials_i = num2cell(1:numTrials)';
        Trials = [Trials; Trials_i];
        
        % get the file name
        filedetails_i = repmat(datalisting(i),numTrials,1);
        filedetails = [filedetails; filedetails_i];
        
        % get the run
        if isempty(strfind(filedetails_i{1},'run'))
            disp('the filename does not contains "run" information')
            return
        end
        run_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'run')+3:strfind(filedetails_i{1},'run')+4)); %run
        runr = [runr; num2cell(repmat(run_thisT,numTrials,1))];
        
        % get the session
        if isempty(strfind(filedetails_i{1},'sess'))
            disp(['--- the filename does not contains "session" information ---'])
            return
        end
        session_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'sess')+4)); %session
        session = [session; num2cell(repmat(session_thisT,numTrials,1))];
        
        % get the std of the prior from the file nm;
        if isempty(strfind(filedetails_i{1},'Pstd'))
            disp(['--- the filename does not contains "Pstd" information ---'])
            return
        end
        Pstd_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2)); %std of the prior
        
        % get the std of the prior from the file nm when prior is uniform;
        if isempty(strfind(filedetails_i{1},'inf'))==0
            Pstd_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2)); %std of the prior
        end
        
        % check if the Prior's std was found in the name
        if isnan(Pstd_thisT)
            disp(['--- the code does not find the value of "Pstd" in the filename ---'])
            return
        end
        Pstd = [Pstd; num2cell(repmat(Pstd_thisT,numTrials,1))];
        
        %get the prior's mean
        if isempty(strfind(filedetails_i{1},'mean'))
            disp(['--- the filename does not contains "Prior mean" information ---'])
            return
        end
        priormean_thisT = str2double(filedetails_i{1}(strfind(filedetails_i{1},'mean')+4:strfind(filedetails_i{1},'mean')+6)); %run
        priormean = [priormean; num2cell(repmat(priormean_thisT,numTrials,1))];
        
        % get the sample directions (group 1)
        % when design is balanced across coherences (same trial number).
        if isfield(datai{1,2}.parameter,'dir')
            sample_dir_i = num2cell(datai{1,2}.parameter.dir)';
            sample_dir   = [sample_dir; sample_dir_i];
            % when design is adjusted according to coherences (different trial nb).
        elseif isfield(datai{1,2}.randVars,'myRandomDir')
            sample_dir_i = num2cell(datai{1,2}.randVars.myRandomDir)';
            sample_dir   = [sample_dir; sample_dir_i];
        end
        % get the sample coherences (group 2)
        % when design is balanced across coherences (same trial number).
        if isfield(datai{1,2}.parameter,'coherence')
            coh_i = num2cell(datai{1,2}.parameter.coherence)';
            coh  = [coh; coh_i];
            % when design is adjusted according to coherences (different trial nb).
        elseif isfield(datai{1,2}.randVars,'myRandomCoh')
            coh_i = num2cell(datai{1,2}.randVars.myRandomCoh)';
            coh  = [coh; coh_i];
        end
    end
    % toc
    % create our databank
    % enter the names of the variable in the columns
    databank.nm = [
        {'filedetails'  },...
        {'run'          },...
        {'session'      },...
        {'Trials'       },...
        {'est_coor'},...
        {'Pstd'         },...
        ('priormean'    ),...
        {'sample_dir'   },...
        {'coh'          },...
        {'est_dir' }];
    
    % store the data in each column
    databanki.data = [
        filedetails(:)  ...
        runr(:)         ...
        session(:)      ...
        Trials(:)       ...
        est_coor(:)     ...
        Pstd(:)         ...
        priormean(:)    ...
        sample_dir(:)   ...
        coh(:)          ...
        est_dirf4Visu];
    
    % store data over subjects
    databank.data = [databank.data; databanki.data];
end

% Data are organized and saved in a file called datbank in the directory
% switch back to the mother directory
cd(pathFold)


% Discard a particular condition from the databank
FAoflvl2deli = strcmp(databank.nm,{FAoflvl2del});
lvl2deli = cell2mat(databank.data(:, FAoflvl2deli)) == lvl2del;
databank.data(lvl2deli,FAoflvl2deli) ={[]};

% backup the databank
save ('datbank','databank');

%% Organize the data (Sort the factors of the analysis)
function [inputLSfitting, FA] = initFAs(FAs, fig)
global databank

% Get est data (e.g., estimated directions)
est_data    = cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));

% Set factors
% FA 1
FA.g1.thisT = cell2mat(databank.data(:,(strcmp(databank.nm, FAs{1}  ))==1));
FA.g1.nm    = FAs{1};
% FA 2
FA.g2.thisT = cell2mat(databank.data(:,(strcmp(databank.nm, FAs{2}  ))==1));
FA.g2.nm    = FAs{2};
% FA 3
FA.g3.thisT = cell2mat(databank.data(:,(strcmp(databank.nm, FAs{3}  ))==1));
FA.g3.nm    = FAs{3};


% Get and order LVLS of each group
% g1 (e.g., priors)
FA.g1.lvlsnm=unique(FA.g1.thisT); %names
FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend'); %order
FA.g1.lvlsnb=numel(FA.g1.lvlsnm);   %number
clear i
for i=1:FA.g1.lvlsnb
    indX.g1.lvli(i) = {find(FA.g1.thisT==FA.g1.lvlsnm(i))};
end
% g2 (e.g., coherences)
FA.g2.lvlsnm=unique(FA.g2.thisT); %names
FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend'); %order
FA.g2.lvlsnb=numel(FA.g2.lvlsnm);   %number
for i=1:FA.g2.lvlsnb
    indX.g2.lvli(i) = {find(FA.g2.thisT==FA.g2.lvlsnm(i))};
end
% g3 (e.g., displayed directions)
FA.g3.lvlsnm=unique(FA.g3.thisT); %names
FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'ascend'); %order
FA.g3.lvlsnb=numel(FA.g3.lvlsnm);   %number
for i=1:FA.g3.lvlsnb
    indX.g3.lvli(i) = {find(FA.g3.thisT==FA.g3.lvlsnm(i))};
end

% Calculate coordinates of average estimated directions for each condition
% organize data in following order: group 1(subplots) - group 2(colors) -
% group 3(x-axis). Each cell contains repetitions of a condition.
% change g1
for k=1:FA.g1.lvlsnb
    %change g2
    for j=1:FA.g2.lvlsnb
        % change g3
        for i=1:FA.g3.lvlsnb
            % indX of single conditions
            indX.g1g2g3(j,i,k) = ...
                {intersect( ...
                intersect ( indX.g1.lvli{k}, indX.g2.lvli{j} ),...
                indX.g3.lvli{i} )};
            % calculate statistics of estimated directions
            est_dir_stat{j,i,k} = statcircular(est_data(indX.g1g2g3{j,i,k},:));
            % extract data
            est_dir.datamean(j,i,k) = {est_dir_stat{j,i,k}.deg.mean};
            % extract information about data
            % g1
            est_dir.taskParams1{j,i,k} = FA.g1.lvlsnm(k);
            % g2
            est_dir.taskParams2{j,i,k} = FA.g2.lvlsnm(j);
            % g3
            est_dir.taskParams3{j,i,k} = FA.g3.lvlsnm(i);
        end
    end
end


% Input for model fitting
inputLSfitting = [est_dir.datamean(:),...
    est_dir.taskParams1(:),...
    est_dir.taskParams2(:),...
    est_dir.taskParams3(:)];

%% Fit the model with the data
function [minSSELLH, freePara2, minSSE, modelPred, udata, taskParams, freeParamsBkp, Rsquared, SSE_bkp, k_bkp] = LSfitting(inputLSfitting, freePara2)

disp('running LSfitting now----------------------------------------------')

% Fit data for each prior separately.
% Parameters search (model fitting by Least Square Optimization (LSO))

% e.g., of the parameters used:
% c=0.12;
% s=2;
% sPr=20;
% uPr=225;
% uLl=155:10:295;
% % one free parameter
% k=s/sPr;

k               = [];
alpha           = [];
beta            = [];
n               = [];
udata           = [];
taskParams1     = [];
taskParams2     = [];
taskParams3     = [];
taskParams      = [];
freeParamsBkp   = [];

% remove "NaN" data
inputLSfitting(cellfun(@(inputLSfitting) any(isnan(inputLSfitting(:,1))),...
    inputLSfitting),:) = [];

% set the condition(s) we want to model
Cond2model.thisT = cell2mat(inputLSfitting(:,3));
Cond2model.lvlnm = unique(Cond2model.thisT);
Cond2model.lvlnm = sort(Cond2model.lvlnm,'descend');
Cond2model.lvlnb = numel(Cond2model.lvlnm);

% options = optimset('MaxFunEvals',1000,'MaxIter',1000);

clear ii
% loop over the conditions
for ii = 1 : Cond2model.lvlnb
    % extract the data
    udata{ii}       = cell2mat(inputLSfitting(Cond2model.thisT == Cond2model.lvlnm(ii) , 1));
    % extract the task parameters
    % parameter 1(e.g., coherence)
    taskParams1{ii} = cell2mat(inputLSfitting(Cond2model.thisT == Cond2model.lvlnm(ii) ,2));
    % parameter 2(e.g., prior's strength)
    taskParams2{ii} = cell2mat(inputLSfitting(Cond2model.thisT == Cond2model.lvlnm(ii) ,3));
    % parameter 3(e.g., displayed direction)
    taskParams3{ii} = cell2mat(inputLSfitting(Cond2model.thisT == Cond2model.lvlnm(ii) ,4));
    % store the parameters
    taskParams {ii}  = [taskParams1{ii}  taskParams2{ii}  taskParams3{ii}];
    
    % Display the parameter 2 of the task
    disp(Cond2model.lvlnm(ii))
    
    
    % Constrained fitting (k cannot be < 0; duration: ?)
    % Set initial values values of the parameters
    
    
    % loop over initial k to avoid local minima
    % k0 = 3;
    k0 = [10 20 40 80];%:1:1000; %exp(-10:0.01:10);
    
%     tic
    clear i
    % Change the initial values of the parameters
    for i = 1 : numel(k0)
%         tic
        % outputs
        SSE         = [];
        freeParams  = [];
        
        % Search the optimal parameters
        [freeParams] = fmincon( @(freeParams) makeSSEI(udata{ii},...
            taskParams{ii},...
            freeParams,...
            freePara2), ...
            k0(i),...
            [],[],...
            [],[],...
            0   ,... 
            +inf,...
            [])%, options);
        
        % Calculate SSE
        SSE = makeSSEI(udata{ii}, taskParams{ii}, freeParams, freePara2);
        % Backup SSE
        SSE_bkp{ii}(i) = SSE;
        % Backup the optimal k for each initial k tested
        k_bkp{ii}(i) = freeParams;
%         toc
    end
  
    
    
    %--------------------------------------------------------------------
    % Draw and store minSSE, model's predictions for minSSE and data
    %--------------------------------------------------------------------
    % value of the lower SSE
    minSSE{ii} = min(SSE_bkp{ii});
    [idxrow{ii}, idxcol{ii}] = min(SSE_bkp{ii});
  
    % parameter for the lower SSE
    k{ii} = k_bkp{ii}(idxcol{ii});
    
    % Record free parameters's info
    freeParamsBkp{1,ii} = Cond2model.lvlnm(ii);
    
    % Compute the model's prediction
    [minSSE{ii}, modelPred{ii}, udata{ii}, freeParamsBkp{2,ii}] = makeSSEI(udata{ii},taskParams{ii},k{ii}, freePara2);
%     [minSSE{ii}, modelPred{ii}, udata{ii}] = makeSSEI(udata{ii},taskParams{ii},k{ii}, freePara2);
    % calculate R^2, the percent of variance in the data explained by the
    % model.
    [Rsquared{ii}] = makeRsquared(udata{ii}, minSSE{ii});
            
%     % Verified the convergence of the optimization
%     figure('color',[1 1 1]);
%     subplot(211); plot(k0, SSE_bkp{ii},'k'); 
%     hold on; 
%     plot(k0(idxcol{ii}), SSE_bkp{ii}(idxcol{ii}), 'O', 'markerfacecolor', 'r', 'markersize', 15)
%     subplot(212); plot(k0, k_bkp{ii},'k'); 
%     hold on; 
%     plot(k0(idxcol{ii}),k_bkp{ii}(idxcol{ii}), 'O', 'markerfacecolor', 'r', 'markersize', 15)    
end

% store output for the optimization of the parameter of the LLH
freePara2 = freePara2;
minSSELLH = sum([minSSE{:}]);

%% Draw the data and the model's predictions
function [fig1] = grapDataAndModel(modelPred, udata, taskParams, FA, fig)
global databank

% Draw group1-lvl1 (e.g., prior std=80)
% Enlarge figure for good quality publication
% fig1.hdle=figure('Position', [0 0 1000 400]); % pixels
fig1.hdle = figure('color',[1 1 1]);
fig1.nm = [fig.nm,'_DataAndModel'];

% Organize the data for plotting
taskParams4plot = [];
udata4plot      = [];
modelPred4plot  = [];
for i = 1 : numel(taskParams)
    % data (e.g., estimated directions)
    udata4plot = [udata4plot; udata{i}];
    % model's predictions
    modelPred4plot  = [modelPred4plot; modelPred{i}];
    % task parameters
    taskParams4plot = [taskParams4plot; taskParams{i}];
end

% % Set factors
% FA 1 (e.g., coherence)
clear i
for i = 1 : FA.g1.lvlsnb
    indX.g1.lvl_i(i) = {find(taskParams4plot(:,1) == FA.g1.lvlsnm(i))};
end
% FA 2 (e.g., priors)
clear i
for i = 1 : FA.g2.lvlsnb
    indX.g2.lvl_i(i) = {find(taskParams4plot(:,2) == FA.g2.lvlsnm(i))};
end
% FA 3 (e.g., displayed directions)
clear i
for i = 1 : FA.g3.lvlsnb
    indX.g3.lvl_i(i) = {find(taskParams4plot(:,3) == FA.g3.lvlsnm(i))};
end


% initialize the parameters of the plot
% set the colors
c = colormap;
FA.g2.color = {[0.5 0 0],...
    c(55,:),...
    [1 0.4 0],...
    [1 0.7 0],...
    c(40,:),...
    c(32,:),...
    c(27,:),...
    c(22,:),...
    c(8 ,:),...
    c(5 ,:),...
    c(1 ,:)}; % group 2(e.g., coherences)

if numel(FA.g2.color) < FA.g2.lvlsnb
    disp(['--- You may want to add ',...
        num2str(FA.g2.lvlsnb - numel(FA.g2.color)),...
        'more colors to the color code ---']);
    return
end
% set axes positions
width=1/(FA.g1.lvlsnb+1);
gap=(1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end
%----------------------------------------------------------------------
% Plot data
%----------------------------------------------------------------------

p12 = [];
p13 = [];

% loop over levels of factor 1 (e.g. coherences)
for k = 1 : FA.g1.lvlsnb
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    % draw background
    %---------------------------------------------------------------------
    hold all
    % draw priors' mean
    p1 = plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],[databank.data{1,7} databank.data{1,7}],...
        'b--',...
        'linewidth',3);
    
    % draw cardinal directions
    p2 = plot([0 30   ], [0   0  ],'color',[0 0 0],'linestyle','--');
    p3 = plot([60 120 ], [90  90 ],'color',[0 0 0],'linestyle','--');
    p4 = plot([150 210], [180 180],'color',[0 0 0],'linestyle','--',...
        'linewidth',3);
    p5 = plot([240 300], [270 270],'color',[0 0 0],'linestyle','--',...
        'linewidth',3);
    
    % Note: add linecode here for drawing an arrow indicating initial position
    % of response line.
    
    % draw the ideal performance line
    p10 = plot([FA.g3.lvlsnm], [FA.g3.lvlsnm'],'k--');
    
    % organize data for plotting
    %---------------------------------------------------------------------
    % change the level of factor 2 (e.g.,prior's strength)
    for j = 1 : FA.g2.lvlsnb
        % change the level of factor 3 (displayed directions)
        % coordinates of single conditions
        indX.g1g2(j,k)={intersect( indX.g1.lvl_i{k},indX.g2.lvl_i{j} ) };
        % store data
        fig1.datamean{j,k}  = udata4plot (indX.g1g2{j,k}); %est_dir{j,i,k}.deg.mean;
        % extract information about data
        % g1
        fig1.dataInfog1{j,k} = taskParams4plot(indX.g1g2{j,k},1);%FA.g1.lvlsnm(k);
        % g2
        fig1.dataInfog2{j,k} = taskParams4plot(indX.g1g2{j,k},2);%FA.g2.lvlsnm(j);
        % g3
        fig1.dataInfog3{j,k} = taskParams4plot(indX.g1g2{j,k},3);%FA.g3.lvlsnm(i);
        
        % store model's predictions
        fig1.modelPred{j,k}  = modelPred4plot(indX.g1g2{j,k});
        
        
        % Normalize data for plotting (because circular data
        % ----------------------------------------------------------------
        % In case the direction is displayed in the 3rd quarter.
        if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
            % "Linearize" the value of the estimated direction
            fig1.datamean{j,k} = 360 + fig1.datamean{j,k};
            % In case the direction is displayed in the 1st quarter.
        elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
            % "Linearize" the value of the estimated direction
            fig1.datamean{j,k} = fig1.datamean{j,k} - 360;
        end
        
        
        % plot data
        %-----------------------------------------------------------------
        p12_ = plot(fig1.dataInfog3{j,k}, fig1.datamean{j,k},'o',... %groups 3 forms x-axis
            'color',FA.g2.color{j},...
            'markerfacecolor',FA.g2.color{j},...
            'markersize',12,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        p12 = [p12 p12_];
        
        % plot model's predictions
        %----------------------------------------------------------------
        p13_ = plot(fig1.dataInfog3{j,k} , [fig1.modelPred{j,k}],'-',... %groups 3 forms x-axis
            'color', FA.g2.color{j},...
            'linewidth', 3,...
            'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        p13 = [p13 p13_];
    end
    leg = legend (p12,'location','northwest');
    set(leg, 'Box', 'off');
    
    % set graph parameters
    % set units of x-axis
    xunit = 1:2:FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'ytick',FA.g3.lvlsnm(xunit),'yticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',12);
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    
    % set labels of axes
    if k==1
        ylabel(ax(k),'est direction','fontsize',14);
    end
    xlabel('sample direction','fontsize',14);
end

% Title
ax(k+1,:) = axes('position',[0 0 1 1],'visible','off');
for k = 1:FA.g1.lvlsnb
    %     text(axs.position(k,1)+axs.position(k,3)/2,...
    %         0.05 + axs.position(k,2)+axs.position(k,4),...
    %         strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
    %         'fontweight','Bold',...
    %         'fontsize',16);
    xpos(k) = axs.position(k,1) + axs.position(k,3)/2;
    ypos(k) = 0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
        'fontweight','Bold',...
        'fontsize',16);
    
end

% Calculate SSE for model II
function [SSE, modelPred, data, freeParams] = makeSSEI(udata, taskParams, freePara, freePara2)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sPr   `: prior's standard deviation (e.g., 20 and inf).
% s      : a constant.
% k:     : a free parameter (k=s/sPr)


% task parameters
%-------------------------------------------------------------------------
% coherence
c   = taskParams(:,1);
% mean of the prior (most likely input)
uPr = 225;
% mean of the likelihood(true input)
uLl = taskParams(:,3);
% free parameters
sp   = freePara;
% sp = 1;


% ------------------------------------------------------------------
% Operations
%-------------------------------------------------------------------
% # in theory, s/c is the variance of the likelihood (i.e., noise) and can't
% be negative
% # in theory sPr is the variance of the prior and can't be negative.
% # What about k = s/sPr ?

% sPr cannot < 0, thus if k<0, it means s<0. If s<0, s/c>0 only if c<0;
% Thus, when k<0, encoding of coherence<0;
% Not sure how that makes sense....

% run Bayesian operation (2nd hypothesis)
% only here changes
%
% The free parameters are 
% #alpha
% #n
% #beta
%-------------------------------------------
% to remove - just for test
alpha = freePara2(1);
% alpha = 1;
beta = freePara2(2);
% beta = 1; 
n = freePara2(3);
% n = 1;


mp = 225;
%------------------------------------------
% Equation when the noise of LLH is a sigmoidal function of the coherence
% sl = alpha./ (c.^n./(c.^n + beta.^n));
% this model is more general. The othe model: sl=s./c corresponds to this
% general model with n=1, alpha=1 and if we add 1;
sl = alpha.*(1 + (beta./c).^n);

% General model
% sl cannot be = inf (defined between ]0:+inf])
% sp cannot be = inf (defined between ]0:+inf])
% uPo = (uLl./sl.^2 + mp.*sp.^2) ./ (1./sl.^2 + 1./sp.^2);
uPo = uLl.*(1./(1 + (sl./sp).^2)) + mp.*(1./(1 + (sp./sl).^2));

%------------------------------------------------
% test - to remove
% uPo=1;
% sp=1;
%------------------------------------------------

% Calculate the SSE between the data and the model's prediction
SSE = sum( (udata - uPo).^2 );

% Store model's informations
modelPred  = uPo;
freeParams = sp;
data       = udata;








% Calculate R^2
function [Rsquared] = makeRsquared(data, SSE)
% Rsquared = 1 - SSE/SST
% SSE (alias residual) is the SSE calculated for the best model's prediction.
% note: R^ (alias coefficient of determination) can be <0 when the model
% does an awful job predicting the data. In this case SSE exceeds that is
% the model fits the data even worse than does a horizontal line.
% Model may not be appropriate or constraints may not be set correctly.
% see http://www.graphpad.com/support/faqid/711/
%
SST = sum((data - mean(data)).^2);
Rsquared = 1 - SSE/SST;

% Calculate the circular statistics of the data
function [data] = statcircular(coord)
% input a vector of cartesian coordinates (coord)
% output :
%   -

% register the coordinates of the input directions
data.coord.all=coord;

% convert from cartesian coordinates to angles (in degree)
data.deg.all = getangle(coord(:,1), coord(:,2));

% calculate the cartesian coordinates of the mean direction est
data.coord.mean = nanmean(coord,1);

% calculate the mean direction est (in degree)
data.deg.mean = getangle(data.coord.mean(:,1),data.coord.mean(:,2));


% calculate the std to the mean direction est (in degree); !!! could be a
% subfunction itself
% Apply the rule of thumb that follows. It seems to work fine intuitively. It would be nice to
% fine a cleaner way to calculate the std.
% initialize the 'sample' and 'mean' variables used to calculate the std
data.num=numel(data.deg.all); % sample size
data.deg.allforstd=data.deg.all;
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

% if the resulting mean direction is between 0 and 180.
if data.deg.mean + 180 <= 360
    % sample each est direction
    for i = 1:data.num
        % if the est direction sampled is >= mean direction + 180
        if data.deg.all(i) >= data.deg.mean + 180
            data.deg.allforstd(i)=data.deg.all(i) - 360;
        end
    end
    % if the resulting mean direction is between 180 and 360.
else
    % sample each est direction
    for i = 1:data.num
        % if the est direction sampled is <= the mean direction - 180
        if data.deg.all(i) <= data.deg.mean - 180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

% now calculate the variance of the est direction.
data.deg.var=nanmean((data.deg.allforstd - data.deg.meanforstd).^2,1);

% and now calculate the std
data.deg.std=sqrt(data.deg.var);

% and now calculate the sem
data.deg.sem=data.deg.std/sqrt(data.num);

% Convert from cartesian coordinates to angles (degree)
function [angle] = getangle(x,y)
% check! to check if the function works fine, write
% e.g., input=180; output=getangle(cos(angle*pi/180),sin(angle*pi/180));
% if the function works input and output should always be the same between
% 0 and 360.

% convert from cartesian coordinates to angle in radians
angle=atan(y./x); % theta=arctan(opposite/adjacent);

% adjust each angle according to his quadrant (in degree)
for i=1:numel(x) % sample each angle
    if x(i)>=0 && y(i)>=0                   %(quadrant 1)
        angle(i) = angle(i)*180/pi;
    elseif x(i)<0                      %(quadrant 2 & 3)
        angle(i) = angle(i)*180/pi + 180;
    elseif x(i)>=0 && y(i)<0               %(quadrant 4)
        angle(i) = angle(i)*180/pi + 360;
    end
end












