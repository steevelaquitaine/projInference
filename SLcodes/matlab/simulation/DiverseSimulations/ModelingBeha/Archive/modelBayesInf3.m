


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


%-----------------------------------------------------------------------------
% Updates
%-----------------------------------------------------------------------------
% 130415 - Normalization of the circular data estimated direction =
% distance (displayed, estimated) + displayed. 
% This manipulation linearizes the data.

% 130419 -  Fit the model to mean and variance of the data pooled together
% I cannot derive the equations for the posterior variance and mean such that I
% get only one free parameter.
% I use the task design (independent manipulation of the prior's variance
% and likelihood's variance) to fit the model with two free parameters:
% prior's variance and likelihood'variance.


%--------------------------------------------------------------------------
% To do
%--------------------------------------------------------------------------
% 130511 - Clean up and improve BI3-5....
% 130516 - add descriptive statistics




function [Sp, modelPred, freePa, Rsquared, udata, sdata, data_for_output, modelPred_for_output, FAbkup, FA, StdParams] = modelBayesInf3(data, FAs, fig, FAoflvl2del, lvl2del )

% OUTPUTS
% ---------------------------------------------------------------------------
modelPred               = []; % Look below: present twice : need to clean up !!!!!!!!!!!!! 
freePa                  = [];
Rsquared                = [];
udata                   = []; % Look below: present twice : need to clean up !!!!!!!!!!!!! 
sdata                   = [];
data_for_output         = []; % Look below: present twice : need to clean up !!!!!!!!!!!!! 
modelPred_for_output    = []; % model prediction restricted to data space
FAbkup                  = []; 
FA                      = [];
Sp                      = [];
StdParams               = [];

% Store data and factors
% ---------------------------------------------------------------------------
% Organize databank
display('Create data matrix')
Makedatabank(data, FAoflvl2del, lvl2del);

% Get factors
display('Get factors')
[inputLSfitting, FA] = initFAs(FAs,fig);


% Modeling
% ---------------------------------------------------------------------------
% % Model BI.1
%--------------------
% % Fit data's mean with width of llh(std) an inverse function of coherence.
% % You are free to set (or not, []) 1 free parameters.
% freePa = [];
% % freePa = [0.1 0.2 0.3 0.4];
% [modelPred, udata, sdata, FAbkup, freePa, Rsquared, data_for_output, modelPred_for_output] = LSfitting(inputLSfitting, freePa);%Fit
% plotModel(modelPred, udata, FAbkup, FA,fig); %Plot model predictions
% Sp = GetPrior1(freePa, FA.g2.lvlsnm); %Get prior's representation

% ----- Try to fit also to std of the data here -----!



% % Model BI.2a
%--------------------
% % Fit data's mean and std with std of llh(std) a hyperbolic function of
% % coherence.
%     % - Bayesian inference
%     % - 4 free parameters for hyperbolic function, std of priors and motor noise
% display('Now fitting the model....')
% % % You are free to set (or not, []) 9 free parameters.
% freePa = [];
% % freePa = [0.1 0.2 0.3 0.4];
% [modelPred, udata, sdata, FAbkup, freePa, Rsquared, data_for_output, modelPred_for_output]  = LSfitting2a(inputLSfitting, freePa );
% display('Now drawing the predictions....')
% plotModel5(modelPred, udata, sdata, FAbkup, FA, fig)
% Sp = GetPrior2a(freePa, FA.g2.lvlsnm); %Get prior's representation



% Model BI.3
%--------------------
%%% Fit data's mean with width of llh(std) as a free parameter and assuming the true priors
% [modelPred, udata, FAbkup, freePa, Rsquared] = LSfitting3(inputLSfitting );
% ----- Try to fit also to std of the data here -----!



% Model BI.4
%--------------------
%%% Fit data's mean and std with width of llh(std) as a free parameter and assuming the true priors
% [modelPred, udata, sdata, FAbkup, freePa, Rsquared] = LSfitting4(inputLSfitting );
% [fig1] = plotModel4(modelPred, udata, sdata, FAbkup, FA,fig);



% % Model BI.5 (~30 min)
%--------------------
%%% Fit data's mean and std with 
    % % - Bayesian inference
    % % - std of llh,  std of priors and motor noise as 8 free parameters
fprintf(' \n %s \n', 'Now fitting the model with the mean and std of the data ...')
%You are free to set (or not, []) your own free parameters.
%freePa = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.01];
freePa=[];
[modelPred, udata, sdata, FAbkup, freePa, Rsquared,data_for_output, modelPred_for_output, StdParams] = LSfitting5(inputLSfitting, freePa);
fprintf(' \n %s \n', 'Now drawing model predictions....')
fig1=plotModel5(modelPred, udata, sdata, FAbkup, FA,fig);
Sp=GetPrior5(freePa, StdParams, FA.g2.lvlsnm, FA.g1.lvlsnm);



%Model BI.6 (~30 min)
%--------------------
%Fit data's mean and std with Bayesian inference
    %std llh,  std priors and motor noise as 8 free parameters
fprintf(' \n %s \n', 'Now fitting the model with the mean and std of the data ...')




%% Create a matrix of cells in which columns are the variables (databank)
function Makedatabank(data, FAoflvl2del,lvl2del)
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

%%Loop over the subfolders to analyse
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
        
        % get the Displayed directions (group 1)
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
            % Find the unique conditions
            indX.g1g2g3(j,i,k) = ...
                {intersect( ...
                intersect ( indX.g1.lvli{k}, indX.g2.lvli{j} ),...
                indX.g3.lvli{i} )};
            
            % Calculate the statistics of the estimated directions
            est_dir_stat{j,i,k} = statcircular(est_data(indX.g1g2g3{j,i,k},:));
            
            % Collect the data
            % Extract the average estimated directions in degree (no normalization)
            est_dir.datamean(j,i,k) = {est_dir_stat{j,i,k}.deg.mean};
            est_dir.datastd(j,i,k) = {est_dir_stat{j,i,k}.deg.std};
            
            % Extract the coordinates of the average estimated directions 
            est_dir.datameanCoord{j,i,k} = est_dir_stat{j,i,k}.coord.mean;
            
            % Collect the information (factors & levels) about the data
            % g1 (e.g., coherences)
            est_dir.FAbkup1{j,i,k} = FA.g1.lvlsnm(k);
            % g2 (e.g., prior)
            est_dir.FAbkup2{j,i,k} = FA.g2.lvlsnm(j);
            % g3 (e.g., displayed directions)
            est_dir.FAbkup3{j,i,k} = FA.g3.lvlsnm(i);
            
            % Collect the displayed directions in cartesian coordinates
            r = 2.5; %the radius of the patch of random dots
            disp.coord{j,i,k} = polar2cartesian(est_dir.FAbkup3{j,i,k},r);
            
            
            % Normalize (linearize) the data for fitting & plotting (because data are circular)
            % ---------------------------------------------------------------------------------
            % Geometric method
            % Calculate the distance between the displayed and the estimated directions.
            est_dir.DistanceDatamean2Displayed{j,i,k} = vectors2signedAngle(est_dir.datameanCoord{j,i,k}, disp.coord{j,i,k});
            % Linearize
            est_dir.datamean_degLinear(j,i,k) = est_dir.FAbkup3{j,i,k} + est_dir.DistanceDatamean2Displayed{j,i,k};
%             % Use the data that have been normalized for the analyses
%             est_dir.datamean(j,i,k) = {est_dir.datamean_degLinear(j,i,k)}; 
            % Use the raw data
            est_dir.datamean(j,i,k) = est_dir.datamean(j,i,k); 
        end
    end
end

% Input for model fitting
inputLSfitting = [est_dir.datamean(:),...%(e.g., mean estimated direction)
    est_dir.datastd(:),...%(e.g., std of estimated direction)
    est_dir.FAbkup1(:),...%(e.g., coherences)
    est_dir.FAbkup2(:),...%(e.g., prior)
    est_dir.FAbkup3(:)]; %(e.g., displayed directions)

% [fig1] = plotModel(modelPred, udata, FAbkup, FA,fig);
%% Fit the model to the data
function [modelPred, udata, sdata, FAbkup, freePa, Rsquared, data_for_output, modelPred_for_output] = LSfitting(inputLSfitting, freePa)
% GOAL
    % Fit mean with one free parameter
    % Least Square Optimization (LSO)

% INPUTS
% "inputLSfitting": matrix of data and factors
    % col1: data
    % col2: factor 1
    % col3: factor 2
    % col4: factor 3
% freePa: vector of free parameters or empty matrix 
    % e.g., freePa = [0.1 0.2 0.3 0.4];
    % e.g., freePa = [];

    
% Get data & factors
%--------------------------------------------------------------------

% Remove "NaN" data
l = inputLSfitting(:,1);
inputLSfitting(cellfun(@(l) any(isnan(l)), l), :) = [];

% Collect data
udata = [inputLSfitting{:,1}]';
sdata = [inputLSfitting{:,2}]';

% Factors
% g1
FA.g1 = [inputLSfitting{:,3}]';
FA.g1lvlnm = unique(FA.g1);
FA.g1lvlnm = sort(FA.g1lvlnm,'descend'); %order

% g2
FA.g2 = [inputLSfitting{:,4}]';
FA.g2lvlnm = unique(FA.g2);
FA.g2lvlnm = sort(FA.g2lvlnm,'descend'); %order

% g3
FA.g3 = [inputLSfitting{:,5}]';
FA.g3lvlnm = unique(FA.g3);
FA.g3lvlnm = sort(FA.g3lvlnm,'ascend'); %order

% Store
FAbkup = [FA.g1 FA.g2 FA.g3];


% If there is no parameters, fit the model to the data
if isempty(freePa)==1
    
    % Initial parameters
    k0_p1 = 0.05:0.05:0.15;
    k0_p2 = 0.05:0.05:0.15;
    k0_p3 = 0.05:0.05:0.15;
    k0_p4 = 0.05:0.05:0.15;
    
    % Fit
    %--------------------------------------------------------------------
    % iteration = 0;
    % Loop over free parameters
    for i = 1 : numel(k0_p1)%e.g., prior 1
        for j = 1 : numel(k0_p2)%e.g., prior 2
            for k = 1 : numel(k0_p3)%e.g., prior 3
                for l = 1 : numel(k0_p4)%e.g., prior 4
                    
                    % Fit
                    freePatmp = fmincon( @(freePatmp) makeSSE1(udata, FAbkup,...
                        freePatmp), ...
                        [k0_p1(i); k0_p2(j); k0_p3(k); k0_p4(l)],...
                        [], [], [], [],...
                        [0 0 0 0],...
                        [+inf +inf +inf +inf], ...
                        []);
                    
                    % Store
                    freePabkp{i,j,k,l} = freePatmp;
                    
                    % Calculate SSE
                    SSE_bkp(i,j,k,l)  = makeSSE1(udata, FAbkup, freePatmp);
                    
                    %     % Check
                    %     fprintf('%6.2f \n',SSE_bkp(i,j,k,l,m,n,e,f))
                    %     iteration = iteration + 1;
                    %     hold all
                    %     plot(iteration, SSE_bkp(i,j,k,l,m,n,e,f),'o', 'markerfacecolor','b',...
                    %         'markersize',13)
                    %     drawnow
                end
            end
        end
    end
    
    % Get the lower SSE
    [minSSE, position] = min(SSE_bkp(:));
    [i, j, k, l] = ind2sub(size(SSE_bkp), position);
    
    % Get the best parameters
    freePa = freePabkp{i, j, k, l};
    
    % Get the model's predictions
    [modelPred, Factors] = makePrediction1(FA, freePa);
    
    % Get the R^2
    % % method 1
    Rsquared = makeRsquared([udata; sdata], minSSE);
    % % method 2
    % [~, modelPred2] = makeSSE5(udata, sdata, FAbkup, freePa);
    % Rsquared = makeRsquared2([udata; sdata], [modelPred2.mean'; modelPred2.std']);
    
    % Store output
    data_for_output = udata;
    modelPred_for_output = modelPred.mean;
    
    % If free parameters are input, only draw predictions
elseif isempty(freePa)==0
    
    % Calculate the SSE
    [~, modelPred00]  = makeSSE1(udata, FAbkup, freePa);
    
    % Calculate the Rsquared
    Rsquared = makeRsquared2(udata, modelPred00.mean');
    
    % Store output
    data_for_output = udata;
    modelPred_for_output = modelPred00.mean';
    
    % Get model predictions on a larger space
    modelPred = makePrediction1(FA, freePa);
    
end

% [fig1] = plotModel(modelPred, udata, FAbkup, FA,fig);
%% Fit the model to the data
function [modelPred, udata, sdata, FAbkup, freePa, Rsquared, data_for_output, modelPred_for_output] = LSfitting2a(inputLSfitting, freePa)
% GOAL
    % Fit a model where the width of the llh is a hyperbolic function of coherence.
        % Busse, L., Ayaz, A., Dhruv, N. T., Katzner, S., Saleem, A. B., Scholvinck, M. L., et al. (2011). 
        % The Detection of Visual Contrast in the Behaving Mouse. Journal of Neuroscience, 31(31), 11351?11361. 
        % doi:10.1523/JNEUROSCI.6689-10.2011
        % Albrecht, D. G., & Hamilton, D. B. (1982). Striate cortex of monkey and cat: contrast response 
        % function. J Neurophysiol.
        % Least Square Optimization (LSO)

% INPUTS
% "inputLSfitting": matrix of data and factors
    % col1: data
    % col2: factor 1
    % col3: factor 2
    % col4: factor 3
% freePa: vector of free parameters or empty matrix 
    % e.g., freePa = [0.1 0.2 0.3 0.4];
    % e.g., freePa = [];

    
% Get data & factors
%--------------------------------------------------------------------

% Remove "NaN" data
l = inputLSfitting(:,1);
inputLSfitting(cellfun(@(l) any(isnan(l)), l), :) = [];

% Collect data
udata = [inputLSfitting{:,1}]';
sdata = [inputLSfitting{:,2}]';

% Factors
% g1
FA.g1 = [inputLSfitting{:,3}]';
FA.g1lvlnm = unique(FA.g1);
FA.g1lvlnm = sort(FA.g1lvlnm,'descend'); %order

% g2
FA.g2 = [inputLSfitting{:,4}]';
FA.g2lvlnm = unique(FA.g2);
FA.g2lvlnm = sort(FA.g2lvlnm,'descend'); %order

% g3
FA.g3 = [inputLSfitting{:,5}]';
FA.g3lvlnm = unique(FA.g3);
FA.g3lvlnm = sort(FA.g3lvlnm,'ascend'); %order

% Store
FAbkup = [FA.g1 FA.g2 FA.g3];

% If there is no parameters, fit the model to the data
if isempty(freePa)==1
    
    % Initial parameters
    R0_0 = 0%:0.05:0.1;
    Rmax_0 = 1%:0.05:0.1;
    n_0 = 1%:0.05:0.1;
    C50_0 = 1%:0.05:0.1;
    sp0_1 = 1:81:164;
    sp0_2 = 1:81:164;
    sp0_3 = 1:81:164;
    sp0_4 = 1:81:164;
    sM0 = 0%:0.05:0.1;

    % Fit
    %--------------------------------------------------------------------
    % iteration = 0;
    % Loop over free parameters
    for i = 1 : numel(R0_0)%e.g.,
        for j = 1 : numel(Rmax_0)%e.g.,
            for k = 1 : numel(n_0)%e.g., motor noise
                for l = 1 : numel(C50_0)%e.g.,
                    for m = 1 : numel(sp0_1)%e.g., prior 1
                        for n = 1 : numel(sp0_2)%e.g., prior 2
                            for o = 1 : numel(sp0_3)%e.g., prior 3
                                for p = 1 : numel(sp0_4)%e.g., prior 4
                                    for q = 1 : numel(sM0)%e.g., motor noise

                                    % Fit
                                    freePatmp = fmincon( @(freePatmp) makeSSE2a(udata, sdata, FAbkup,...
                                        freePatmp), ...
                                        [R0_0(i); Rmax_0(j); n_0(k); C50_0(l); sp0_1(m); sp0_2(n); sp0_3(o); sp0_4(p); sM0(q)],...
                                        [], [], [], [],...
                                        [0 0 0 0 0 0 0 0 0],...
                                        [+inf +inf +inf +inf +inf +inf +inf +inf +inf], ...
                                        []);
                                    
                                    % Store
                                    freePabkp{i,j,k,l,m,n,o,p,q} = freePatmp;
                                    
                                    % Calculate SSE
                                    SSE_bkp(i,j,k,l,m,n,o,p,q)  = makeSSE2a(udata, sdata, FAbkup, freePatmp);
                                    
                                    %     % Check
                                    %     fprintf('%6.2f \n',SSE_bkp(i,j,k,l,m,n,e,f))
                                    %     iteration = iteration + 1;
                                    %     hold all
                                    %     plot(iteration, SSE_bkp(i,j,k,l,m,n,e,f),'o', 'markerfacecolor','b',...
                                    %         'markersize',13)
                                    %     drawnow
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Get the lower SSE
    [minSSE, position] = min(SSE_bkp(:));
    [i,j,k,l,m,n,o,p,q] = ind2sub(size(SSE_bkp), position);
    
    % Get the best parameters
    freePa = freePabkp{i,j,k,l,m,n,o,p,q};
    
    % Get the model's predictions
    [modelPred, ~] = makePrediction2a(FA, freePa);
    
    % Get the R^2
    % % method 1
    Rsquared = makeRsquared([udata; sdata], minSSE);
    % % method 2
    % [~, modelPred2] = makeSSE5(udata, sdata, FAbkup, freePa);
    % Rsquared = makeRsquared2([udata; sdata], [modelPred2.mean'; modelPred2.std']);
    
    % Store output
    data_for_output = udata;
    modelPred_for_output = modelPred.mean;
    
    % If free parameters are input, only draw predictions
elseif isempty(freePa)==0
    
    % Calculate the SSE
    [~, modelPred00]  = makeSSE1(udata, FAbkup, freePa);
    
    % Calculate the Rsquared
    Rsquared = makeRsquared2(udata, modelPred00.mean');
    
    % Store output
    data_for_output = udata;
    modelPred_for_output = modelPred00.mean';
    
    % Get model predictions on a larger space
    modelPred = makePrediction2a(FA, freePa);
    
end


%% Fit the model to the data
function [modelPred, udata, FAbkup, freePa, Rsquared, fig1] = LSfitting3(inputLSfitting )

% Fit data for each prior separately.
% Parameters search (model fitting by Least Square Optimization (LSO))

% e.g., of the parameters used:
% c=0.12;
% s=2;
% sPr=20;
% uPr=225;
% uLl=155:10:295;
% % two free parameters:
%  sl, sp


% Remove "NaN" data
l=inputLSfitting(:,1);
inputLSfitting(cellfun(@(l) any(isnan(l)),...
    l),:) = [];


% Fit sl for each coherence
% for data(coh == 1), set the priors as true and fit sl.

% Collect the data for coherence ii
FAtofit.labels = cell2mat(inputLSfitting(:,3));%(e.g., coh)
FAtofit.lvlnm = sort(unique(FAtofit.labels),'descend');

% Initial parameters
k0 = 0.05:0.05:20;


freePaBkp   = [];
% Loop over coherences
for ii = 1:numel(FAtofit.lvlnm);

    %--------------------------------------------------------------------
    % Select variables
    %--------------------------------------------------------------------
    % Data
    udata{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii), 1));
    
    % Conditions
    % g1
    FA.g1{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii) ,3));
    % g2
    FA.g2{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii) ,4));
    % g3
    FA.g3{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii) ,5));
    % Store
    FAbkup{ii} = [FA.g1{ii}  FA.g2{ii}  FA.g3{ii}];
   
    %--------------------------------------------------------------------
    % Fit data
    %--------------------------------------------------------------------
    % Loop over the initial parameters
    for i = 1 : numel(k0)
        
        % Fitting
        freePatmp = fmincon( @(freePatmp) makeSSE3(udata{ii}, FAbkup{ii},...
            freePatmp), ...
            k0(i),...
            [], [], [], [], 0, +inf, []);
        
        % Store free parameters
        freePaBkup{1,ii}(i) = freePatmp;
        
        % Calculate the SSE
        SSE_bkp{ii}(i) = makeSSE3(udata{ii}, FAbkup{ii}, freePaBkup{1,ii}(i));

    end
    
    %--------------------------------------------------------------------
    % Draw SSE, predictions and data
    %--------------------------------------------------------------------
    % Get the lower SSE
    minSSE{ii} = min(SSE_bkp{ii});
    [~, idxcol] = min(SSE_bkp{ii});
  
    % Get the parameter for the lower SSE and the factor
    freePa{ii,1} = freePaBkup{ii}(idxcol);
    freePa{ii,2} = FAtofit.lvlnm(ii);
    
    % Compute the model's best predictions
    [~, modelPred{ii}, ~] = makeSSE3(udata{ii}, FAbkup{ii}, freePa{ii,1});
    
    % Calculate R^2, the percent of variance in the data explained by the
    % model.
    Rsquared{ii} = makeRsquared(udata{ii}, minSSE{ii});        
    
    
    % Check the convergence of the optimization
    figure('color', [1 1 1]);
    subplot(211); 
    title('SSE','fontsize',20)
    hold all
    plot(k0, SSE_bkp{ii},'k'); 
    plot(k0(idxcol), minSSE{ii}, 'O', 'markerfacecolor', 'r', 'markersize', 10)
    subplot(212);
    title('sl','fontsize',20)
    hold all
    plot(k0, freePaBkup{1,ii},'k'); 
    plot(k0(idxcol),freePa{ii,1}, 'O', 'markerfacecolor', 'r', 'markersize', 10)    
    
end

% [fig1] = plotModel(modelPred, udata, FAbkup, FA,fig);

%% Fit the model to the data
function [modelPred, udata, sdata, FAbkup, freePa, Rsquared, fig1] = LSfitting4(inputLSfitting )
global fig
% Fit mean and variance with std(llh) as a free parameter and assuming the true priors
% Parameters search (model fitting by Least Square Optimization (LSO))

% e.g., of the parameters used:
% c=0.12;
% s=2;
% sPr=20;
% uPr=225;
% uLl=155:10:295;
% % two free parameters:
%  sl, sp


% Remove "NaN" data
l=inputLSfitting(:,1);
inputLSfitting(cellfun(@(l) any(isnan(l)),...
    l),:) = [];

% Fit sl for each coherence
% for data(coh == 1), set the priors as true and fit sl.

% Collect the data for coherence ii
FAtofit.labels = cell2mat(inputLSfitting(:,3));%(e.g., coh)
FAtofit.lvlnm = sort(unique(FAtofit.labels),'descend');

% Initial parameters
k0 = 0.05%:0.05:20;


freePaBkp   = [];
% Loop over coherences
for ii = 1:numel(FAtofit.lvlnm);

    %--------------------------------------------------------------------
    % Select variables
    %--------------------------------------------------------------------
    % Data
    udata{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii), 1));
    sdata{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii), 2));

    % Conditions
    % g1
    FA.g1{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii) ,3));%e.g.,coh
    % g2
    FA.g2{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii) ,4));
    % g3
    FA.g3{ii} = cell2mat(inputLSfitting(FAtofit.labels == FAtofit.lvlnm(ii) ,5));
    % Store
    FAbkup{ii} = [FA.g1{ii}  FA.g2{ii}  FA.g3{ii}];
   
    %--------------------------------------------------------------------
    % Fit data
    %--------------------------------------------------------------------
    % Loop over the initial parameters
    for i = 1 : numel(k0)
        
        % Fitting
        freePatmp = fmincon( @(freePatmp) makeSSE4(udata{ii}, sdata{ii}, FAbkup{ii},...
            freePatmp), ...
            k0(i),...
            [], [], [], [], 0, +inf, []);
        
        % Store free parameters
        freePaBkup{1,ii}(i) = freePatmp;
        
        % Calculate the SSE
        SSE_bkp{ii}(i) = makeSSE4(udata{ii}, sdata{ii}, FAbkup{ii}, freePaBkup{1,ii}(i));

    end
    
    %--------------------------------------------------------------------
    % Draw SSE, predictions and data
    %--------------------------------------------------------------------
    % Get the lower SSE
    minSSE{ii} = min(SSE_bkp{ii});
    [~, idxcol] = min(SSE_bkp{ii});
  
    % Get the parameter for the lower SSE and the factor
    freePa{ii,1} = freePaBkup{ii}(idxcol);
    freePa{ii,2} = FAtofit.lvlnm(ii);
    
    % Compute the model's best predictions
    [~, modelPred{ii}, ~] = makeSSE4(udata{ii}, sdata{ii}, FAbkup{ii}, freePa{ii,1});
    
    % Calculate R^2, the percent of variance in the data explained by the
    % model.
    Rsquared{ii} = makeRsquared([udata{ii};sdata{ii}], minSSE{ii});        
   
    % Check the convergence of the optimization
    figure('color', [1 1 1]);
    subplot(211); 
    title('SSE','fontsize',20)
    hold all
    plot(k0, SSE_bkp{ii},'k'); 
    plot(k0(idxcol), minSSE{ii}, 'O', 'markerfacecolor', 'r', 'markersize', 10)
    subplot(212);
    title('sl','fontsize',20)
    hold all
    plot(k0, freePaBkup{1,ii},'k'); 
    plot(k0(idxcol),freePa{ii,1}, 'O', 'markerfacecolor', 'r', 'markersize', 10)    
    xlabel('initial values')
end
% [fig1] = plotModel4(modelPred, udata, sdata, FAbkup, FA,fig);


%% Fit the model to the data
function [modelPred, udata, sdata, FAbkup, freePa, Rsquared, data_for_output, modelPred_for_output, StdParams] = LSfitting5(inputLSfitting, freePa)
% Fit mean and variance with std(llh) as a free parameter and assuming the true priors
% Parameters search (model fitting by Least Square Optimization (LSO))

% Remove "NaN" data
l = inputLSfitting(:,1);
inputLSfitting(cellfun(@(l) any(isnan(l)), l), :) = [];

%--------------------------------------------------------------------
% Select variables
%--------------------------------------------------------------------
% Collect data
udata = [inputLSfitting{:,1}]';
sdata = [inputLSfitting{:,2}]';

% Conditions
% g1
FA.g1 = [inputLSfitting{:,3}]';
FA.g1lvlnm = unique(FA.g1);
FA.g1lvlnm = sort(FA.g1lvlnm,'descend'); %order

% g2
FA.g2 = [inputLSfitting{:,4}]';
FA.g2lvlnm = unique(FA.g2);
FA.g2lvlnm = sort(FA.g2lvlnm,'descend'); %order

% g3
FA.g3 = [inputLSfitting{:,5}]';
FA.g3lvlnm = unique(FA.g3);
FA.g3lvlnm = sort(FA.g3lvlnm,'ascend'); %order

% Store
FAbkup = [FA.g1 FA.g2 FA.g3];


% If there is no parameters, fit the model to the data
if isempty(freePa) == 1
    
    % Initial parameters(Elapsed time is 2114.179696 seconds)
    sl1_0 = 1:81:164;%std of likelihood: c=0.24
    sl2_0 = 1:81:164;
    sl3_0 = 1:81:164;
    sp1_0 = 1:81:164;%std of prior: 80 deg
    sp2_0 = 1:81:164;
    sp3_0 = 1:81:164;
    sp4_0 = 1:81:164;
    sM_0  = 0;%1:81:164;%motor noise
    
    %--------------------------------------------------------------------
    % Fitting
    %--------------------------------------------------------------------
    % iteration = 0;
    % Loop over free parameters
    for i = 1 : numel(sl1_0)
        for j = 1 : numel(sl2_0)
            for k = 1 : numel(sl3_0)
                for l = 1 : numel(sp1_0)
                    for m = 1 : numel(sp2_0)
                        for n = 1 : numel(sp3_0)
                            for e = 1 : numel(sp4_0)
                                for f = 1 : numel(sM_0)
                                    
                                    %                                     % Fit
                                    %                                     [freePatmp] = fmincon( @(freePatmp) makeSSE5(udata, sdata, FAbkup,...
                                    %                                         freePatmp), ...
                                    %                                         [sl1_0(i); sl2_0(j); sl3_0(k); sp1_0(l); sp2_0(m); sp3_0(n); sp4_0(e); sM_0(f)],...
                                    %                                         [], [], [], [],...
                                    %                                         [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001],...
                                    %                                         [+inf +inf +inf +inf +inf +inf +inf +inf], []);
                                    %
                                    %Fit
                                    [freePatmp,~,~,~,~,~,Hessian] = fmincon( @(freePatmp) makeSSE5(udata, sdata, FAbkup,...
                                        freePatmp), ...
                                        [sl1_0(i); sl2_0(j); sl3_0(k); sp1_0(l); sp2_0(m); sp3_0(n); sp4_0(e); sM_0(f)],...
                                        [], [], [], [],...
                                        [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001],...
                                        [+inf +inf +inf +inf +inf +inf +inf +inf], []);
                                                                        
                                    % Store
                                    freePabkp{i,j,k,l,m,n,e,f} = freePatmp;
                                    
                                    % Calculate the SSE
                                    SSE_bkp(i,j,k,l,m,n,e,f)  = makeSSE5(udata, sdata, FAbkup, freePatmp);
                                    
%                                     % Check
%                                     fprintf('%6.2f \n',SSE_bkp(i,j,k,l,m,n,e,f))
%                                     iteration = iteration + 1;
%                                     hold all
%                                     plot(iteration, SSE_bkp(i,j,k,l,m,n,e,f),'o', 'markerfacecolor','b',...
%                                         'markersize',13)
%                                     drawnow
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    % Get the lower SSE
    [minSSE, position] = min(SSE_bkp(:));
    [i,j,k,l,m,n,e,f] = ind2sub(size(SSE_bkp),position);
    
    % Get the best parameters
    freePa = freePabkp{i,j,k,l,m,n,e,f};
    
    % Get extended model's predictions
    [modelPred, ~] = makePrediction(FA, freePa);
    
    % Calculate data-restricted model predictions
    [~, modelPred00]  = makeSSE5(udata, sdata, FAbkup, freePa);
    
    % Get the R^2
    % % method 1
    Rsquared = makeRsquared([udata;sdata], minSSE);
    % % method 2
    % [~, modelPred2] = makeSSE5(udata, sdata, FAbkup, freePa);
    % Rsquared = makeRsquared2([udata; sdata], [modelPred2.mean'; modelPred2.std']);
    
    % Store output
    data_for_output = [udata; sdata];
    modelPred_for_output = [modelPred00.mean'; modelPred00.std'];
    
    % Get standard deviation of model parameters
    StdParams = GetStdparams(Hessian, data_for_output, modelPred_for_output, freePa);
    
% If free parameters are input, only draw predictions
elseif isempty(freePa)==0
    
    % Calculate the SSE
    [~, modelPred00]  = makeSSE5(udata, sdata, FAbkup, freePa);
    
    % Calculate the Rsquared
    Rsquared = makeRsquared2([udata; sdata], [modelPred00.mean'; modelPred00.std']);
    
    % Store output
    data_for_output = [udata; sdata];
    modelPred_for_output = [modelPred00.mean'; modelPred00.std'];
    
    % Get the model's predictions generalized
    [modelPred, ~] = makePrediction(FA, freePa);
    
    % Getting standard deviation of model parameters is not meaningful here
    stdParams = []; 
end



%% Draw the data and the model's predictions
function [fig1] = plotModel(modelPred, udata, FAbkup, FA, fig)

%%% GOAL: 
% Draw data sorted based on factors 1, 2 & 3
% A - Plots: data against factor 3 & factor 2
% B - Subplots: data against factor 1

%%% INPUT:
    %   Factors: experimental factors
    %     udata: data
    % modelPred: model predictions

global databank


% Draw group1-lvl1 (e.g., prior std=80)
% Enlarge figure for good quality publication
% fig1.hdle=figure('Position', [0 0 1000 400]); % pixels
fig1.hdle = figure('color',[1 1 1]);
fig1.nm = [fig.nm,'_DataAndModel'];

% % Set factors
% FA 1 (e.g., coherence)
clear i
for i = 1 : FA.g1.lvlsnb
    indX.g1.lvl_i(i) = {find(FAbkup(:,1) == FA.g1.lvlsnm(i))};
end
% FA 2 (e.g., priors)
clear i
for i = 1 : FA.g2.lvlsnb
    indX.g2.lvl_i(i) = {find(FAbkup(:,2) == FA.g2.lvlsnm(i))};
end
% FA 3 (e.g., displayed directions)
clear i
for i = 1 : FA.g3.lvlsnb
    indX.g3.lvl_i(i) = {find(FAbkup(:,3) == FA.g3.lvlsnm(i))};
end

% initialize the parameters of the plot
% set the colors
c = colormap;
FA.g2.color = {[0.5 0 0],...
    c(55,:),...
    [1 0.4 0],...
    [0.75 0.75 0],...
    c(40,:),...
    c(32,:),...model
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
% Plot data - mean
%----------------------------------------------------------------------
% space axes
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% loop over levels of factor 1 (e.g. coherences)
for k = 1 : FA.g1.lvlsnb 
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    % draw background
    %---------------------------------------------------------------------
    hold all
    % draw priors' mean
    p1 = plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],[databank.data{1,7} databank.data{1,7}],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    % Note: add linecode here for drawing an arrow indicating initial position
    % of response line.
    
    % draw the ideal performance line
    p10 = plot(FA.g3.lvlsnm, FA.g3.lvlsnm', 'k:',...
        'linewidth',2,...
        'Displayname','Ideal predictions');
    
    % organize data for plotting
    %---------------------------------------------------------------------
    % change the level of factor 2 (e.g., prior's strength)
    for j = 1 : FA.g2.lvlsnb
        % change the level of factor 3 (displayed directions)
        % coordinates of single conditions
        indX.g1g2(j,k) = {intersect( indX.g1.lvl_i{k},indX.g2.lvl_i{j} ) };
        % store data
        fig1.datamean{j,k} = udata(indX.g1g2{j,k}); %est_dir{j,i,k}.deg.mean;

        % extract information about data
        % g1
        fig1.dataInfog1{j,k} = FAbkup(indX.g1g2{j,k},1);%FA.g1.lvlsnm(k);
        % g2
        fig1.dataInfog2{j,k} = FAbkup(indX.g1g2{j,k},2);%FA.g2.lvlsnm(j);
        % g3
        fig1.dataInfog3{j,k} = FAbkup(indX.g1g2{j,k},3);%FA.g3.lvlsnm(i);

               
%         % Normalize data for plotting (because circular data) -already
%         done earlier in the code. should be suppressed
%         % ----------------------------------------------------------------
%         % In case the direction is displayed in the 3rd quarter.
%         if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = 360 + fig1.datamean{j,k};
%             % In case the direction is displayed in the 1st quarter.
%         elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = fig1.datamean{j,k} - 360;
%         end

       
        % plot data
        %-----------------------------------------------------------------
        p12_ = plot(fig1.dataInfog3{j,k}, fig1.datamean{j,k},'o',... %groups 3 forms x-axis
            'color',FA.g2.color{j},...
            'markerfacecolor',FA.g2.color{j},...
            'markersize',15,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        %p12 = [p12 p12_];
        
        
        % plot model's predictions
        %----------------------------------------------------------------
        % mean
        p13a = plot(FA.g3.lvlsnm, modelPred.mean(:,j, k),'-',... %groups 3 forms x-axis
            'color', FA.g2.color{j} - [0.2 0 0],...
            'linewidth', 2,...
            'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
    end

    % set the graph parameters
    % set the unit steps of the x axis
    xunit = 1:6:FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'ytick',FA.g3.lvlsnm(xunit),'yticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    % Set the limits of the axes
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);  
    % set the labels of the axes
    if k==1
        ylabel(ax(k),'Estimated directions','fontsize',20);
    end
    xlabel('Displayed directions','fontsize',20);
end
% % Legend the graph
% leg = legend (p12,'location','BestOutside');
% set(leg, 'Box', 'off');

% Title
ax(k+1,:) = axes('position',[0 0 0.97 0.75],'visible','off');
for k = 1 : FA.g1.lvlsnb
    %     text(axs.position(k,1)+axs.position(k,3)/2,...
    %         0.05 + axs.position(k,2)+axs.position(k,4),...
    %         strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
    %         'fontweight','Bold',...
    %         'fontsize',20);
    xpos(k) = axs.position(k,1) + axs.position(k,3)/2;
    ypos(k) = 0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
        'fontweight','Bold',...
        'fontsize',20);
    
end

%% Draw the data and the model's predictions
function [fig1] = plotModel4(modelPred, udata, sdata, FAbkup, FA, fig)
global databank

% Draw group1-lvl1 (e.g., prior std=80)
% Enlarge figure for good quality publication
% fig1.hdle=figure('Position', [0 0 1000 400]); % pixels
fig1.hdle = figure('color',[1 1 1]);
fig1.nm = [fig.nm,'_DataAndModel'];

% Organize the data for plotting
FAbkup4plot = [];
udata4plot      = [];
sdata4plot      = [];
modelPred4plotmean  = [];
modelPred4plotstd  = [];

for i = 1 : numel(FAbkup)
    % mean of the data (e.g., estimated directions)
    udata4plot = [udata4plot; udata{i}];
    % sdt of the data (e.g., std estimated directions)
    sdata4plot = [sdata4plot; sdata{i}];
    % model's predictions mean and std
    modelPred4plotmean = [modelPred4plotmean; modelPred{i}.mean];
    modelPred4plotstd  = [modelPred4plotstd; modelPred{i}.std];
    % task parameters
    FAbkup4plot = [FAbkup4plot; FAbkup{i}];
end

% % Set factors
% FA 1 (e.g., coherence)
clear i
for i = 1 : FA.g1.lvlsnb
    indX.g1.lvl_i(i) = {find(FAbkup4plot(:,1) == FA.g1.lvlsnm(i))};
end
% FA 2 (e.g., priors)
clear i
for i = 1 : FA.g2.lvlsnb
    indX.g2.lvl_i(i) = {find(FAbkup4plot(:,2) == FA.g2.lvlsnm(i))};
end
% FA 3 (e.g., displayed directions)
clear i
for i = 1 : FA.g3.lvlsnb
    indX.g3.lvl_i(i) = {find(FAbkup4plot(:,3) == FA.g3.lvlsnm(i))};
end


% initialize the parameters of the plot
% set the colors
c = colormap;
FA.g2.color = {[0.5 0 0],...
    c(55,:),...
    [1 0.4 0],...
    [0.75 0.75 0],...
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
% Plot data - mean
%----------------------------------------------------------------------
% space axes
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% loop over levels of factor 1 (e.g. coherences)
for k = 1 : FA.g1.lvlsnb
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    % draw background
    %---------------------------------------------------------------------
    hold all
    % draw priors' mean
    p1 = plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],[databank.data{1,7} databank.data{1,7}],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    % Note: add linecode here for drawing an arrow indicating initial position
    % of response line.
    
    % draw the ideal performance line
    p10 = plot([FA.g3.lvlsnm], [FA.g3.lvlsnm'],'k:',...
        'linewidth',1.00005,...
        'Displayname','Ideal predictions');
    
    % organize data for plotting
    %---------------------------------------------------------------------
    % change the level of factor 2 (e.g.,prior's strength)
    for j = 1 : FA.g2.lvlsnb
        % change the level of factor 3 (displayed directions)
        % coordinates of single conditions
        indX.g1g2(j,k)={intersect( indX.g1.lvl_i{k},indX.g2.lvl_i{j} ) };
        % store data
        fig1.datamean{j,k} = udata4plot(indX.g1g2{j,k}); %est_dir{j,i,k}.deg.mean;

        % extract information about data
        % g1
        fig1.dataInfog1{j,k} = FAbkup4plot(indX.g1g2{j,k},1);%FA.g1.lvlsnm(k);
        % g2
        fig1.dataInfog2{j,k} = FAbkup4plot(indX.g1g2{j,k},2);%FA.g2.lvlsnm(j);
        % g3
        fig1.dataInfog3{j,k} = FAbkup4plot(indX.g1g2{j,k},3);%FA.g3.lvlsnm(i);
        
        % store model's predictions
        fig1.modelPred.mean{j,k}  = modelPred4plotmean(indX.g1g2{j,k});
        fig1.modelPred.std{j,k}  = modelPred4plotstd(indX.g1g2{j,k});

               
%         % Normalize data for plotting (because circular data) -already
%         done earlier in the code. should be suppressed
%         % ----------------------------------------------------------------
%         % In case the direction is displayed in the 3rd quarter.
%         if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = 360 + fig1.datamean{j,k};
%             % In case the direction is displayed in the 1st quarter.
%         elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = fig1.datamean{j,k} - 360;
%         end

       
        % plot data
        %-----------------------------------------------------------------
        p12_ = plot(fig1.dataInfog3{j,k}, fig1.datamean{j,k},'o',... %groups 3 forms x-axis
            'color',FA.g2.color{j},...
            'markerfacecolor',FA.g2.color{j},...
            'markersize',15,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        %p12 = [p12 p12_];
        
        
        % plot model's predictions
        %----------------------------------------------------------------
        % mean
        p13a = plot(fig1.dataInfog3{j,k} , [fig1.modelPred.mean{j,k}],'-',... %groups 3 forms x-axis
            'color', FA.g2.color{j} - [0.2 0 0],...
            'linewidth', 1.0005,...
            'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
    end

    % set the graph parameters
    % set the unit steps of the x axis
    xunit = 1:6:FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'ytick',FA.g3.lvlsnm(xunit),'yticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    % Set the limits of the axes
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);  
    % set the labels of the axes
    if k==1
        ylabel(ax(k),'Estimated directions','fontsize',20);
    end
    xlabel('Displayed directions','fontsize',20);
end
% % Legend the graph
% leg = legend (p12,'location','BestOutside');
% set(leg, 'Box', 'off');

% Title
ax(k+1,:) = axes('position',[0 0 0.97 0.75],'visible','off');
for k = 1:FA.g1.lvlsnb
    %     text(axs.position(k,1)+axs.position(k,3)/2,...
    %         0.05 + axs.position(k,2)+axs.position(k,4),...
    %         strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
    %         'fontweight','Bold',...
    %         'fontsize',20);
    xpos(k) = axs.position(k,1) + axs.position(k,3)/2;
    ypos(k) = 0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
        'fontweight','Bold',...
        'fontsize',20);
    
end

  


%----------------------------------------------------------------------
% Plot data - std
%----------------------------------------------------------------------

fig2.hdle = figure('color','w');
fig2.nm = [fig.nm, '_std'];

% set the axes positions
width = 1/(FA.g1.lvlsnb+1);

% space axes
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% loop over levels of factor 1 (e.g. coherences)
for k = 1 : FA.g1.lvlsnb
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    % draw background
    %---------------------------------------------------------------------
    hold all
    % draw priors' mean
    p1 = plot([databank.data{1,7} databank.data{1,7}],[0 150],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    % organize data for plotting
    %---------------------------------------------------------------------
    % change the level of factor 2 (e.g.,prior's strength)
    for j = 1 : FA.g2.lvlsnb
        % change the level of factor 3 (displayed directions)
        % coordinates of single conditions
        indX.g1g2(j,k)={intersect( indX.g1.lvl_i{k},indX.g2.lvl_i{j} ) };
        
        % store data
        fig1.datastd{j,k} = sdata4plot(indX.g1g2{j,k}); 
        
        % extract information about data
        % g1
        fig1.dataInfog1{j,k} = FAbkup4plot(indX.g1g2{j,k},1);%FA.g1.lvlsnm(k);
        % g2
        fig1.dataInfog2{j,k} = FAbkup4plot(indX.g1g2{j,k},2);%FA.g2.lvlsnm(j);
        % g3
        fig1.dataInfog3{j,k} = FAbkup4plot(indX.g1g2{j,k},3);%FA.g3.lvlsnm(i);
        
        % Store model's predictions
        fig1.modelPred.std{j,k}  = modelPred4plotstd(indX.g1g2{j,k});

        
%         % Normalize data for plotting (because circular data) -already
%         done earlier in the code. should be suppressed
%         % ----------------------------------------------------------------
%         % In case the direction is displayed in the 3rd quarter.
%         if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = 360 + fig1.datamean{j,k};
%             % In case the direction is displayed in the 1st quarter.
%         elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = fig1.datamean{j,k} - 360;
%         end

       
        % plot data
        %-----------------------------------------------------------------
        p12_ = plot(fig1.dataInfog3{j,k}, fig1.datastd{j,k},'o',... %groups 3 forms x-axis
            'color',FA.g2.color{j},...
            'markerfacecolor',FA.g2.color{j},...
            'markersize',15,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
                
        % plot model's predictions
        %----------------------------------------------------------------
        % mean
        p13a = plot(fig1.dataInfog3{j,k} , [fig1.modelPred.std{j,k}],'-',... %groups 3 forms x-axis
            'color', FA.g2.color{j} - [0.2 0 0],...
            'linewidth', 1.0005,...
            'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
    end

    % set the graph parameters
    % set the unit steps of the x axis
    xunit = 1:6:FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    % Set the limits of the axes
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([0 150]);  
    % set the labels of the axes
    if k==1
        ylabel(ax(k),'Std of estimated directions','fontsize',20);
    end
    xlabel('Displayed directions','fontsize',20);
end
% % Legend the graph
% leg = legend (p12,'location','BestOutside');
% set(leg, 'Box', 'off');
    
% Title
ax(k+1,:) = axes('position',[0 0 0.97 0.75],'visible','off');
for k = 1:FA.g1.lvlsnb
    %     text(axs.position(k,1)+axs.position(k,3)/2,...
    %         0.05 + axs.position(k,2)+axs.position(k,4),...
    %         strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
    %         'fontweight','Bold',...
    %         'fontsize',20);
    xpos(k) = axs.position(k,1) + axs.position(k,3)/2;
    ypos(k) = 0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
        'fontweight','Bold',...
        'fontsize',20);
    
end


% backup figures informations
figs = {fig1, fig2};

%% Draw the data and the model's predictions
function [fig1] = plotModel5(modelPred, udata, sdata, FAbkup, FA, fig)
global databank

% Draw group1-lvl1 (e.g., prior std=80)
% Enlarge figure for good quality publication
% fig1.hdle=figure('Position', [0 0 1000 400]); % pixels
fig1.hdle = figure('color',[1 1 1]);
fig1.nm = [fig.nm,'_DataAndModel'];


% ------!!!!!! Seriously check the correspondence of udata, sdata and
% FAbkup !!!!!!!!!!!!!!!!!!!!!!!!



% ------------------------------------------------------------------
% Set the factors
% ------------------------------------------------------------------
% FA 1 (e.g., coherence)
clear i
for i = 1 : FA.g1.lvlsnb
    indX.g1.lvl_i(i) = {find(FAbkup(:,1) == FA.g1.lvlsnm(i))};
end
% FA 2 (e.g., priors)
clear i
for i = 1 : FA.g2.lvlsnb
    indX.g2.lvl_i(i) = {find(FAbkup(:,2) == FA.g2.lvlsnm(i))};
end
% FA 3 (e.g., displayed directions)
clear i
for i = 1 : FA.g3.lvlsnb
    indX.g3.lvl_i(i) = {find(FAbkup(:,3) == FA.g3.lvlsnm(i))};
end


%----------------------------------------------------------------------
% Plot 
%----------------------------------------------------------------------
% Set the graph's parameters
% colors
c = colormap;
FA.g2.color = {[0.5 0 0],...
    c(55,:),...
    [1 0.4 0],...
    [0.75 0.75 0],...
    c(40,:),...
    c(32,:),...
    c(27,:),...
    c(22,:),...
    c(8 ,:),...
    c(5 ,:),...
    c(1 ,:)}; % group 2(e.g., coherences)

% if not enough colors
if numel(FA.g2.color) < FA.g2.lvlsnb
    disp(['--- You may want to add ',...
        num2str(FA.g2.lvlsnb - numel(FA.g2.color)),...
        'more colors to the color code ---']);
    return
end

% axes' positions
width=1/(FA.g1.lvlsnb+1);
gap=(1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% space axes
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end


% Mean

% loop over levels of factor 1 (e.g. coherences)
for k = 1 : FA.g1.lvlsnb
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    % draw background
    %---------------------------------------------------------------------
    hold all
    % draw priors' mean
    p1 = plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],[databank.data{1,7} databank.data{1,7}],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    % Draw the predictions of ideal performances
    p10 = plot(FA.g3.lvlsnm, FA.g3.lvlsnm','k:',...
        'linewidth',1,...
        'Displayname','Ideal predictions');
    
    % organize data for plotting
    %---------------------------------------------------------------------
    % loop over factor 2's levels (e.g.,prior's strength)
    for j = 1 : FA.g2.lvlsnb
        
        % Loop over factor 3's levels (displayed directions)
        % get the positions of each condition
        indX.g1g2(j,k) = {intersect( indX.g1.lvl_i{k},indX.g2.lvl_i{j} ) };
        % Store data
        fig1.datamean{j,k} = udata(indX.g1g2{j,k}); %est_dir{j,i,k}.deg.mean;

        % extract information about data
        % g1
        fig1.dataInfog1{j,k} = FAbkup(indX.g1g2{j,k},1);%FA.g1.lvlsnm(k);
        % g2
        fig1.dataInfog2{j,k} = FAbkup(indX.g1g2{j,k},2);%FA.g2.lvlsnm(j);
        % g3
        fig1.dataInfog3{j,k} = FAbkup(indX.g1g2{j,k},3);%FA.g3.lvlsnm(i);
              
%         % Normalize data for plotting (because circular data) -already
%         done earlier in the code. should be suppressed
%         % ----------------------------------------------------------------
%         % In case the direction is displayed in the 3rd quarter.
%         if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = 360 + fig1.datamean{j,k};
%             % In case the direction is displayed in the 1st quarter.
%         elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = fig1.datamean{j,k} - 360;
%         end

       
        % plot data
        %-----------------------------------------------------------------
        p12_ = plot(fig1.dataInfog3{j,k}, fig1.datamean{j,k} ,'o',... %groups 3 forms x-axis
            'color',FA.g2.color{j},...
            'MarkerEdgeColor','w',...
            'markerfacecolor',FA.g2.color{j},...
            'markersize',15,...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));        
        
        % plot model's predictions
        %----------------------------------------------------------------
        % mean
        p13a = plot(FA.g3.lvlsnm , modelPred.mean(:,j,k),'-',... %groups 3 forms x-axis
            'color', FA.g2.color{j} - [0.2 0 0],...
            'linewidth', 2,...
            'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        
        % Store data and predictions plotted here....
        
    end

    % set the graph's parameters
    % set the unit steps of the x axis
    xunit = 1: 6: FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'ytick',FA.g3.lvlsnm(xunit),'yticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    % Set the limits of the axes
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);  
    % set the labels of the axes
    if k==1
        ylabel(ax(k),'Estimated directions','fontsize',20);
    end
    xlabel('Displayed directions','fontsize',20);
end
% % Legend the graph
% leg = legend (p12,'location','BestOutside');
% set(leg, 'Box', 'off');

% Title
ax(k+1,:) = axes('position',[0 0 0.97 0.75],'visible','off');
for k = 1:FA.g1.lvlsnb
    xpos(k) = axs.position(k,1) + axs.position(k,3)/2;
    ypos(k) = 0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
        'fontweight','Bold',...
        'fontsize',20); 
end

  


%----------------------------------------------------------------------
% Plot data - std
%----------------------------------------------------------------------

fig2.hdle = figure('color','w');
fig2.nm = [fig.nm, '_std'];

% set the axes positions
width = 1/(FA.g1.lvlsnb+1);

% space axes
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% loop over levels of factor 1 (e.g. coherences)
for k = 1 : FA.g1.lvlsnb
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    % draw background
    %---------------------------------------------------------------------
    hold all
    % draw priors' mean
    p1 = plot([databank.data{1,7} databank.data{1,7}],[0 150],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    % organize data for plotting
    %---------------------------------------------------------------------
    % Loop over factor 2's levels
    for j = 1 : FA.g2.lvlsnb
        
        % Loop over factor 3's levels
        % Get the position of each conditions
        indX.g1g2(j,k)={intersect( indX.g1.lvl_i{k},indX.g2.lvl_i{j} ) };
        
        % Collect the data
        fig1.datastd{j,k} = sdata(indX.g1g2{j,k}); 
        
        % Extract the factors/levels label
        % g1
        fig1.dataInfog1{j,k} = FAbkup(indX.g1g2{j,k},1);%FA.g1.lvlsnm(k);
        % g2
        fig1.dataInfog2{j,k} = FAbkup(indX.g1g2{j,k},2);%FA.g2.lvlsnm(j);
        % g3
        fig1.dataInfog3{j,k} = FAbkup(indX.g1g2{j,k},3);%FA.g3.lvlsnm(i);
        
%         % Normalize data for plotting (because circular data) -already
%         done earlier in the code. should be suppressed
%         % ----------------------------------------------------------------
%         % In case the direction is displayed in the 3rd quarter.
%         if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = 360 + fig1.datamean{j,k};
%             % In case the direction is displayed in the 1st quarter.
%         elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
%             % "Linearize" the value of the estimated direction
%             fig1.datamean{j,k} = fig1.datamean{j,k} - 360;
%         end

       
        % plot the data
        %-----------------------------------------------------------------
        p12_ = plot(fig1.dataInfog3{j,k}, fig1.datastd{j,k},'o',... %groups 3 forms x-axis
            'color',FA.g2.color{j},...
            'markerfacecolor',FA.g2.color{j},...
            'markersize',15,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
                
        % plot the model's predictions
        %----------------------------------------------------------------
        % mean
        p13a = plot(FA.g3.lvlsnm , modelPred.std(:,j,k)','-',... %groups 3 forms x-axis
            'color', FA.g2.color{j} - [0.2 0 0],...
            'linewidth', 2,...
            'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
    end

    % set the graph parameters
    % set the unit steps of the x axis
    xunit = 1:6:FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    % Set the limits of the axes
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([0 150]);  
    % set the labels of the axes
    if k==1
        ylabel(ax(k),'Std of estimated directions','fontsize',20);
    end
    xlabel('Displayed directions','fontsize',20);
end
% % Legend the graph
% leg = legend (p12,'location','BestOutside');
% set(leg, 'Box', 'off');
    
% Title
ax(k+1,:) = axes('position',[0 0 0.97 0.75],'visible','off');
for k = 1:FA.g1.lvlsnb
    xpos(k) = axs.position(k,1) + axs.position(k,3)/2;
    ypos(k) = 0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(FA.g1.nm,': ',num2str(FA.g1.lvlsnm(k))),...
        'fontweight','Bold',...
        'fontsize',20);
end


% backup figures informations
figs = {fig1, fig2};






%% Nested functions
% Model fitting
% Calculate SSE for model I
function [SSE, modelPred, freePa] = makeSSE1(udata, FAbkup, freePa)
% % ------------------------------------------------------------------
% % Representations
% %-------------------------------------------------------------------
% % c      : coherence (e.g., 12,35,100%).
% % uLl    : displayed direction.
% % uPr    : prior's mean, i.e., most likely direction for gaussians.
% % uPo    : posterior's mean, i.e., estimated direction.
% % sPr   `: prior's standard deviation (e.g., 20 and inf).
% % s      : a constant.
% % k:     : a free parameter (k=s/sPr)
% 
% % task parameters
% % Coherence
% c   = FAbkup(:,1);
% % Mean of the prior (most likely input)
% uPr = 225;
% % Mean of the likelihood(true input)
% uLl = FAbkup(:,3);
% 
% 
% % Free parameter
% k   = freePara;
% 
% % ------------------------------------------------------------------
% % Operations
% %-------------------------------------------------------------------
% % # in theory, s/c is the variance of the likelihood (i.e., noise) and can't
% % be negative
% % # in theory sPr is the variance of the prior and can't be negative.
% % # What about k = s/sPr ?
% 
% % sPr cannot < 0, thus if k<0, it means s<0. If s<0, s/c>0 only if c<0;
% % Thus, when k<0, encoding of coherence<0;
% % Not sure how that makes sense....
% 
% % Run Bayesian inference (1st hypothesis, mean of the data)
% uPo = uLl.*( 1./(1+(k./c).^2) ) + uPr.*( 1./(1+(c./k).^2) );
% 
% % Calculate the SSE between the data and the model's prediction
% SSE = sum( (udata - uPo).^2 );
% 
% % Store the model's informations
% modelPred.mean  = uPo;
% freePa = k;
% data       = udata;
% 


% Task parameters
% Coherence
c    = FAbkup(:,1); %coherence
spe  = FAbkup(:,2); %std of the prior
uLl  = FAbkup(:,3); %Mean of the likelihood
uPr  = 225; %mean of the prior

% Debuging
% fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
%     [sl1;sl2;sl3;sp1;sp2;sp3;sp4])

% Assign the free parameters to each condition
for i = 1 : numel(udata)
    
    % Set free parameters
    if spe(i) == 80
        k(i) = freePa(1);
    end
    
    if spe(i) == 40
        k(i) = freePa(2);
    end
    
    if spe(i) == 20
        k(i) = freePa(3);
    end
    
    if spe(i) == 10
        k(i) = freePa(4);
    end
    
    % Run Bayesian inference
    % https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
    % mean of the posterior
    uPo(i) = ( 1./(1+(k(i)./c(i)).^2) ).*uLl(i)  +  ( 1./(1+(c(i)./k(i)).^2) ).*uPr;

end


% Calculate the SSE
SSE = sum( (udata - uPo').^2 );

% fprintf('%6.2f  %12.8f\n',[udata;sdata])
% fprintf('%6.2f  %12.8f\n',[uPo;sPo]);
%         
% Store the model's informations
modelPred.mean = uPo;

% Calculate SSE for model 3
function [SSE, modelPred, freePa] = makeSSE2a(udata, sdata, FAbkup, freePa)



% Representations
%-------------------------------------------------------------------
% Task parameters
c    = FAbkup(:,1); %coherence
spe  = FAbkup(:,2); %std of the prior
uPr = 225;% Mean of the prior
uLl = FAbkup(:,3);% Mean of the likelihood

% Free parameters
%llh
R0  = freePa(1);
Rmax = freePa(2);
n   = freePa(3);
C50 = freePa(4);
%prior
sp1 = freePa(5);%std of prior 1
sp2 = freePa(6);  
sp3 = freePa(7);  
sp4 = freePa(8); 
%motor
sM  = freePa(9);%motor noise

%------------------------------------------
% Equation when the noise of LLH is an hyperbolic function of the coherence
    % Busse, L., Ayaz, A., Dhruv, N. T., Katzner, S., Saleem, A. B., Scholvinck, M. L., 
    % et al. (2011). The Detection of Visual Contrast in the Behaving Mouse. Journal 
    % of Neuroscience, 31(31), 11351?11361. doi:10.1523/JNEUROSCI.6689-10.2011

    % sl = 1./(R0 + Rmax.*((c^n)./(C50 + (c^n))))
    % This model is more general than the previous model.

% Assign the free parameters to each condition
for i = 1 : numel(udata)
    
    % Model width of the likelihood
    sl(i) = 1./(R0 + Rmax.*((c(i).^n)./(C50 + (c(i).^n))));
    
    % Set prior conditions
    if spe(i) == 80
        sp(i) = sp1;
    end
    
    if spe(i) == 40
        sp(i) = sp2;
    end
    
    if spe(i) == 20
        sp(i) = sp3;
    end
    
    if spe(i) == 10
        sp(i) = sp4;
    end
    
    % Bayesian inference
    % mean of the posterior
    uPo(i) = uLl(i).*(1./(1 + (sl(i)./sp(i)).^2)) + (1./(1 + (sp(i)./sl(i)).^2)).*uPr;
    % std of the posterior
    sPo(i) = sqrt(1./((1./sl(i).^2) + (1./sp(i).^2)));
    % std of the estimate
    sEs(i) = sPo(i) + sM;

end

% Calculate the SSE
SSE = sum( ([udata; sdata] - [uPo'; sEs']).^2 );

% fprintf('%6.2f  %12.8f\n',[udata;sdata])
% fprintf('%6.2f  %12.8f\n',[uPo;sPo]);
%         
% Store the model's informations
modelPred.mean = uPo;
modelPred.std  = sEs;
freePa = [R0;Rmax;n;C50;sp1;sp2;sp3;sp4;sM];

% Calculate SSE for model 3
function [SSE, modelPred, freePa] = makeSSE3(udata, FAbkup, freePara)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sl     : std of the likelihood.
% sPr   `: std of the prior (e.g., 20 and inf).
% s      : a constant.
% k:     : a free parameter (k=s/sPr)

% Task parameters
% Coherence
c    = FAbkup(:,1); %coherence
sp   = FAbkup(:,2); %std of the prior
% Mean of the prior
uPr = 225;
% Mean of the likelihood
uLl = FAbkup(:,3);
% Free parameter
sl   = freePara;

% Run Bayesian inference to predict the mean of the data
% https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
uPo = (1./( 1+(sl./sp).^2 )).*uLl + (1./( 1+(sp./sl).^2 )).*uPr;
% sPo = 1./((1./sl.^2) + (1./sp.^2));

% Calculate the SSE
SSE = sum( (udata - uPo).^2 );

% Store the model's informations
modelPred.mean  = uPo;
freePa = sl;
data       = udata;

% Calculate SSE for model 4
function [SSE, modelPred, freePa] = makeSSE4(udata, sdata, FAbkup, freePara)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sl     : std of the likelihood.
% sPr   `: std of the prior (e.g., 20 and inf).
% s      : a constant.
% k:     : a free parameter (k=s/sPr)

% Task parameters
% Coherence
c    = FAbkup(:,1); %coherence
sp   = FAbkup(:,2); %std of the prior
% Mean of the prior
uPr = 225;
% Mean of the likelihood
uLl = FAbkup(:,3);
% Free parameter
sl   = freePara;

% Run Bayesian inference to predict the mean of the data
% https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
% mean
uPo = (1./( 1+(sl./sp).^2 )).*uLl + (1./( 1+(sp./sl).^2 )).*uPr;
% std
sPo = sqrt(1./((1./sl.^2) + (1./sp.^2)));

% Calculate the SSE
SSE = sum( ([udata;sdata] - [uPo;sPo]).^2 );

% Store the model's informations
modelPred.mean = uPo;
modelPred.std  = sPo;
freePa = sl;

% Calculate SSE for model 5
function [SSE, modelPred, freePa] = makeSSE5(udata, sdata, FAbkup, freePa)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sl     : std of the likelihood.
% sPr   `: std of the prior (e.g., 20 and inf).
% s      : a constant.

% udata  :  
% sdata  : 

% Task parameters
c    = FAbkup(:,1); %coherence
spe  = FAbkup(:,2); %std of the prior
uLl  = FAbkup(:,3); %Mean of the likelihood
uPr  = 225; %mean of the prior

% Set the free parameters
sl1 = freePa(1);%std of likelihood 1
sl2 = freePa(2);
sl3 = freePa(3); 
sp1 = freePa(4);%std of prior 1
sp2 = freePa(5);  
sp3 = freePa(6);  
sp4 = freePa(7); 
sM  = freePa(8);%motor noise

% Debuging
% fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
%     [sl1;sl2;sl3;sp1;sp2;sp3;sp4])

% Assign the free parameters to each condition
for i = 1 : numel(udata)
    
    % set the conditions
    if c(i) == 0.24
        sl(i) = sl1;
    end
    
    if c(i) == 0.12
        sl(i) = sl2;
    end
    
    if c(i) == 0.06
        sl(i) = sl3;
    end
    
    if spe(i) == 80
        sp(i) = sp1;
    end
    
    if spe(i) == 40
        sp(i) = sp2;
    end
    
    if spe(i) == 20
        sp(i) = sp3;
    end
    
    if spe(i) == 10
        sp(i) = sp4;
    end
    
    % Run Bayesian inference
    % https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
    % mean of the posterior
    uPo(i) = (1./( 1+(sl(i)./sp(i)).^2 )).*uLl(i) + (1./( 1+(sp(i)./sl(i)).^2 )).*uPr;
    % std of the posterior
    sPo(i) = sqrt(1./((1./sl(i).^2) + (1./sp(i).^2)));
    % std of the estimate
    sEs(i) = sPo(i) + sM;
end


% Calculate the SSE
SSE = sum( ([udata;sdata] - [uPo'; sEs']).^2 );

% fprintf('%6.2f  %12.8f\n',[udata;sdata])
% fprintf('%6.2f  %12.8f\n',[uPo;sPo]);
%         
% Store the model's informations
modelPred.mean = uPo;
modelPred.std  = sEs;
freePa = [sl1;sl2;sl3;sp1;sp2;sp3;sp4;sM];

% Calculate the prediction for model 5
function [modelPred, FA] = makePrediction1(FA, BestfreePa)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sl     : std of the likelihood.
% sPr   `: std of the prior (e.g., 20 and inf).
% s      : a constant.

% Model fixed parameters
c    = FA.g1lvlnm; %coherence
uPr  = 225; %mean of the prior
uLl = FA.g3lvlnm;

% Model free parameters
k = BestfreePa;

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

% Assign the free parameters to each condition
% loop over 
for i = 1 : numel(FA.g1lvlnm) %e.g.,coherence 
    for j = 1 : numel(FA.g2lvlnm) %e.g., prior
        for l = 1 : numel(FA.g3lvlnm) %e.g., displayed directions
            
            % Run Bayesian inference (1st hypothesis, mean of the data)
            uPo(l,j,i) =  ( 1./(1+(k(j)./c(i)).^2) ).*uLl(l)  +  uPr.*( 1./(1+(c(i)./k(j)).^2) );
            
        end
    end
end

% Store the model's informations
modelPred.mean = uPo;

% Calculate the prediction for model 2a
function [modelPred, FA] = makePrediction2a(FA, BestfreePa)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sl     : std of the likelihood.
% sPr   `: std of the prior (e.g., 20 and inf).
% s      : a constant.

% udata  :  
% sdata  : 

% Representations
%-------------------------------------------------------------------
% Task parameters
c   = FA.g1lvlnm; %coherence
uPr = 225;% Mean of the prior
uLl = FA.g3lvlnm;% Mean of the likelihood

% Free parameters
%llh
R0  = BestfreePa(1);
Rmax = BestfreePa(2);
n   = BestfreePa(3);
C50 = BestfreePa(4);
%prior
sp(1) = BestfreePa(5);%std of prior 1
sp(2) = BestfreePa(6);  
sp(3) = BestfreePa(7);  
sp(4) = BestfreePa(8); 
%motor
sM  = BestfreePa(9);%motor noise

% Debuging
% fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
%     [sl1;sl2;sl3;sp1;sp2;sp3;sp4])

% Assign the free parameters to each condition
% loop over 
for i = 1 : numel(FA.g1lvlnm)%e.g.,coherence 
    for j = 1 : numel(FA.g2lvlnm) %e.g., prior
        for k = 1 : numel(FA.g3lvlnm)
            
            % Run Bayesian inference
            % https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
            
            % width of the likelihood
            sl(i) = 1./(R0 + Rmax.*((c(i).^n)./(C50 + (c(i).^n))));

            % mean of the posterior
            uPo(k,j,i) = (1./( 1+(sl(i)./sp(j)).^2 )).*uLl(k) + (1./( 1+(sp(j)./sl(i)).^2 )).*uPr;
            % std of the posterior
            sPo(k,j,i) = sqrt(1./((1./sl(i).^2) + (1./sp(j).^2)));
            % std of the estimate
            sEs(k,j,i) = sPo(k,j,i) + sM;
            
            
            % Store conditions
            FA.g1stored(k,j,i) = FA.g1lvlnm(i);
            FA.g2stored(k,j,i) = FA.g2lvlnm(j);
            FA.g3stored(k,j,i) = FA.g3lvlnm(k);
        end
    end
end

% Store the model's informations
modelPred.mean = uPo;
modelPred.std  = sEs;

% Calculate the prediction for model 5
function [modelPred, FA] = makePrediction(FA, BestfreePa)
% ------------------------------------------------------------------
% Representations
%-------------------------------------------------------------------
% c      : coherence (e.g., 12,35,100%).
% uLl    : displayed direction.
% uPr    : prior's mean, i.e., most likely direction for gaussians.
% uPo    : posterior's mean, i.e., estimated direction.
% sl     : std of the likelihood.
% sPr   `: std of the prior (e.g., 20 and inf).
% s      : a constant.

% udata  :  
% sdata  : 

% Task parameters
uPr  = 225; %mean of the prior

% Input
uLl = FA.g3lvlnm;

% Set the free parameters
sl(1) = BestfreePa(1);%std of likelihood 1
sl(2) = BestfreePa(2);
sl(3) = BestfreePa(3); 
sp(1) = BestfreePa(4);%std of prior 1
sp(2) = BestfreePa(5);  
sp(3) = BestfreePa(6);  
sp(4) = BestfreePa(7); 
sM    = BestfreePa(8);%motor noise

% Debuging
% fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
%     [sl1;sl2;sl3;sp1;sp2;sp3;sp4])

% Assign the free parameters to each condition
% loop over 
for i = 1 : numel(FA.g1lvlnm)%e.g.,coherence 
    for j = 1 : numel(FA.g2lvlnm) %e.g., prior
        for k = 1 : numel(FA.g3lvlnm)
            
            % Run Bayesian inference
            % https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
            % mean of the posterior
            uPo(k,j,i) = (1./( 1+(sl(i)./sp(j)).^2 )).*uLl(k) + (1./( 1+(sp(j)./sl(i)).^2 )).*uPr;
            % std of the posterior
            sPo(k,j,i) = sqrt(1./((1./sl(i).^2) + (1./sp(j).^2)));
            % std of the estimate
            sEs(k,j,i) = sPo(k,j,i) + sM;
            
            
            % Store conditions
            FA.g1stored(k,j,i) = FA.g1lvlnm(i);
            FA.g2stored(k,j,i) = FA.g2lvlnm(j);
            FA.g3stored(k,j,i) = FA.g3lvlnm(k);
        end
    end
end

% Store the model's informations
modelPred.mean = uPo;
modelPred.std  = sEs;

% Calculate R^2 (the variance explained by the model)
function Rsquared = makeRsquared(data, SSE)
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

% Calculate R^2 (the variance explained by the model)
function Rsquared = makeRsquared2(data, modelPred)
% Rsquared = 1 - SSE/SST
% SSE (alias residual) is the SSE calculated for the best model's prediction.
% note: R^ (alias coefficient of determination) can be <0 when the model
% does an awful job predicting the data. In this case SSE exceeds that is
% the model fits the data even worse than does a horizontal line.
% Model may not be appropriate or constraints may not be set correctly.
% see http://www.graphpad.com/support/faqid/711/
%
R = corr(data, modelPred);
Rsquared = R^2;

% Get standard deviation of model parameters
function StdParams = GetStdparams(Hessian, data, modelPred, freePa)
% note: complexe number Std means parameters have negative variance

% Calculate parameters' covariance matrix
% d.covar = inv(jacobian'*jacobian); 
d.covar = inv(Hessian);  %a matrix

% Calculate noise variance
residual = [data - modelPred]';
noiseVariance = (residual*residual')/(numel(data) -  numel(freePa)); % a scalar

% Calculate std of model parameters
StdParams = diag(sqrt(noiseVariance * d.covar))';

% Get the prior's representation
function Sp = GetPrior1(freeP, factor)

% INPUTS
%freeP: In descending order
%factor: correspond to the freeP

% Convert cell input to matrix
if iscell(freeP) == 1
    freeP = [freeP{:}];
end

% Remove cases of infinite prior from the analysis
if freeP(:,1) == +inf 
    freeP(1) = [];
end

% Get the priors' true and estimated values
Sp.exp.raw = factor;
Sp.estimated.raw =  freeP;

% Get normalized priors' width (std of weaker prior/ std of prior_i)
for i = 1 : numel(freeP )
    % Experimental
    Sp.exp.normalized(i) = Sp.exp.raw(1)/Sp.exp.raw(i); %Sw/Ss
    
    % Estimated 
    Sp.estimated.normalized(i) = Sp.estimated.raw(i)/Sp.estimated.raw(1); %Ks/Kw = Sw/Ss
end


% Show summary results
% Variables
fprintf('\n\n\n\n\n %15s %15s \n',...
    'Std.exp.norm',...
    'Std.est.norm')
% Values
for i = 1: numel(Sp.exp.normalized)
    fprintf('\n %15i %15f \n',...
        [Sp.exp.normalized(i)'; Sp.estimated.normalized(i)']);
end

% Draw subjects' data (i.e., the strength of the prior as perceived by subjects)
figure('color',[1 1 1]); 
hold all
title ('Subject representations of prior width')
xlabel('{\sigma}_w_e_a_k_e_r _p_r_i_o_r/{\sigma}_p_r_i_o_r', 'fontsize',20)
ylabel({'K_p_r_i_o_r/K_w_e_a_k_e_r _p_r_i_o_r',...
    'i.e., {\sigma}_w_e_a_k_e_r _p_r_i_o_r/{\sigma}_p_r_i_o_r'}, 'fontsize',20)

plot(Sp.exp.normalized, Sp.estimated.normalized,'-ko',...
    'markerfacecolor','k',... 
    'markersize', 15,...
    'displayname','subject' );
plot(Sp.exp.normalized, Sp.exp.normalized,'k:',...
    'linewidth',6,'displayname', 'True ratio');

ylim([1 10])
xlim([0 max(Sp.exp.normalized)+2])
lg=legend('location','Northwest');
box(lg,'off')
axis square
set(gca,'fontsize',20)

% Get the prior's representation
function Sp = GetPrior2a(freeP, factor)

% INPUTS
    %freeP: In descending order
    %factor: factors associated to the freeP

% Convert cell input to matrix
if iscell(freeP) == 1
    freeP = [freeP{:}];
end

% Remove cases of infinite prior from the analysis
if freeP(:,1) == +inf 
    freeP(1) = [];
end

% Representation of the prior
%---------------------------------
% Get the priors' true and estimated values
Sp.exp.raw        = factor;
Sp.estimated.raw  =  freeP(5:8);

% Show summary results
% Variables
fprintf('\n\n\n\n\n %15s %15s \n',...
    'Std.exp.norm',...
    'Std.est.norm')
% Values
for i = 1: numel(Sp.exp.raw)
    fprintf('\n %15i %15f \n',...
        [Sp.exp.raw(i)'; Sp.estimated.raw(i)']);
end

% Draw subjects' data (i.e., the strength of the prior as perceived by subjects)
figure('color',[1 1 1]); 
hold all
title ('Subject representations of prior width')
xlabel('{\sigma}_p_r_i_o_r', 'fontsize',20)
ylabel('{\sigma}_p_r_i_o_r', 'fontsize',20)

plot(Sp.exp.raw, Sp.estimated.raw,'-ko',...
    'markerfacecolor','k',... 
    'markersize', 15,...
    'displayname','subject' );
plot(Sp.exp.raw, Sp.exp.raw,'k:',...
    'linewidth',6,'displayname', 'True ratio');

% ylim([1 10])
xlim([0 max(Sp.exp.raw)+2])
lg=legend('location','Northwest');
box(lg,'off')
axis square
set(gca,'fontsize',20)
%%% -------  TO DO : add a linear fit to subjects' data !!!!!!
%%% -------  plot also the width of the likelihood !!!!!!!!


% Representation of the likelihood
%---------------------------------
% Get the width of the likelihood
R0      = freeP(1);
Rmax    = freeP(2);
n       = freeP(3);
C50     = freeP(4);

% Set factor 1, e.g., coherence
c = [0.06 0.12 0.24];
% loop over levels of factor 1
for i = 1 : numel(c)
    sl(i) = 1./(R0 + Rmax.*((c(i).^n)./(C50 + (c(i).^n))));
end
figure('color','w'); 
hold all
plot(c, sl,'k-o', 'markerfacecolor','k','markersize',15,'linewidth',2)
ylabel('Std of the likelihood','fontsize', 20)
xlabel('Coherence', 'fontsize', 20)
set(gca, 'fontsize', 20)

% Get the prior's representation
function Sp = GetPrior5(freeP, stdParams, factor, factor2)

% INPUTS
    %freeP: In descending order
    %factor: factors associated to the freeP

    
% Representation of the prior
%---------------------------------
% Get the priors' true and estimated values
Sp.exp.raw        = factor;
Sp.estimated.raw  = freeP(4:7);
Sp.std            = stdParams(4:7);

% Draw subjects' data (i.e., the strength of the prior as perceived by subjects)
figure('color',[1 1 1]); 
hold all
title ('Subject representations of prior width','fontsize',20)
xlabel('Experimental {\sigma}_p_r_i_o_r', 'fontsize',20)
ylabel('Estimated {\sigma}_p_r_i_o_r', 'fontsize',20)
myerrorbar(Sp.exp.raw, Sp.estimated.raw,'yError', Sp.std,...
    'Symbol=o',...
    'Markersize=30',...
    'Color=[0 0 0]');
plot(Sp.exp.raw, Sp.exp.raw, 'k:', 'linewidth', 1)
linefit(Sp.exp.raw, Sp.estimated.raw, 'k')
set(gca,'xtick',[0:10:90],'xticklabel',[0:10:90], 'fontsize', 20)

% Representation of the likelihood
%---------------------------------
% Get the llh' true and estimated values
Sl.estimated.raw  = freeP(1:3);
Sl.std            = stdParams(1:3);

% Draw subjects' data (i.e., the strength of the prior as perceived by subjects)
figure('color', [1 1 1]); 
hold all
title ('Subject representations of llh width','fontsize',20)
xlabel('Coherence', 'fontsize',20)
ylabel('Estimated {\sigma}_l_l_h', 'fontsize',20)
myerrorbar(factor2, Sl.estimated.raw, 'yError', Sl.std,'Color=[0 0 0]');
xlim([0 1])
ylim([0 Sl.estimated.raw(end)+Sl.std(end)])
set(gca,'xtick',[0.06:0.06:1],'xticklabel', [0.06:0.06:1], 'fontsize', 20)


%% Nested functions
% Basic functions
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
    % sample each Estimated direction
    for i = 1:data.num
        % if the Estimated direction sampled is >= mean direction + 180
        if data.deg.all(i) >= data.deg.mean + 180
            data.deg.allforstd(i)=data.deg.all(i) - 360;
        end
    end
    % if the resulting mean direction is between 180 and 360.
else
    % sample each Estimated direction
    for i = 1:data.num
        % if the Estimated direction sampled is <= the mean direction - 180
        if data.deg.all(i) <= data.deg.mean - 180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

% now calculate the variance of the Estimated direction.
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

% Convert from polar to cartesian coordinates
function [coord] = polar2cartesian(theta,r)
% theta is an angle in degree
% r is the radius of the unit circle
% Coord are in visual angle
% Record angle in degree
theta2.deg = theta;
% Convert from degree to radian
theta2.rad = theta2.deg*pi/180;
% Calculate visual angles coordinates
x = r*cos(theta2.rad);
y = r*sin(theta2.rad);
coord = [x y];

% Calculate the angle formed by two vectors (distanceRsquared
function [angle] = vectors2signedAngle(v1,v2)
% Inputs are two vectors' coordinates
% xV1 = v1(1)
% yV1 = v1(2)
% xV2 = v2(1)
% yV2 = v2(2)

%e.g., v1.x=0; v1.y=1; v2.x=1;v2.y=0;
%angle = - (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x, v1.x*v2.x+v1.y*v2.y)
% gives 90 degrees.

%v1.x=1; v1.y=0; v2.x=0;v2.y=1;
%angle = - (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x, v1.x*v2.x+v1.y*v2.y)
% gives - 90 degrees.

% load data
xV1 = v1(1);
yV1 = v1(2);
xV2 = v2(1);
yV2 = v2(2);

% Calculate the angle in degree separating the two vectors
angle = - (180/pi) * atan2(xV1*yV2 - yV1*xV2, xV1*xV2+yV1*yV2);

% Fit the data with a linear model
function linefit(x,y,linecolor)
% find non missing data
datahere=~isnan(y);
y=y(datahere);

% compute x and y for the linear fit
xfit=x;
P = polyfit(x(datahere),y,1);
yfit = polyval(P,x);
plot(xfit,yfit,'-',...
    'linewidth',0.5,...
    'color',linecolor);

% Plot errorbars
function retval = myerrorbar(x,y,varargin)

% Plot errorbar
% myerrorbar.m
%
%      usage: myerrorbar(x,y,varargin)
%         by: justin gardner
%       date: 06/24/07
%    purpose: draw plots with error bars
%       e.g.: y error bars
%             myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10));
%             x error bars
%             myerrorbar(1:10,rand(1,10),'xError',0.5*rand(1,10));
%             x and yerror bars
%             myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10),'xError',0.5*rand(1,10));
%             different lower and upper bounds
%             myerrorbar(1:10,rand(1,10),'yLow',2*rand(1,10),'yHigh',0.5*rand(1,10));
%    options: Symbol = symbol to use, default 'o-'
%             Color = symbol color, default 'k'
%             MarkerFaceColor = symbol face color, defaults to Color
%             MarkerEdgeColor = symbol edge color, defaults to Color
%             MarkerSize = symbol size, default 8
%             myerrorbar(1:10,rand(1,10),'yError',rand(1,10)/2,'Symbol=s-','Color=[1 0.5 0]');
%             tee = draw tees or not, default 0
%             yTeelen = length of tee on y error, default to 1/10 of x spacing
%             xTeelen = length of tee on x error, default to 1/10 of y spacing

% check arguments
if nargin < 2
  help myerrorbar
  return
end
 
% check for old style usage
if (nargout == 1) || ((length(varargin) >= 1) && isnumeric(varargin{1}))
  retval = myerrorbarold(x,y,varargin);
  return
end

% get arguments
getArgs(varargin);

% no passed in x
if ieNotDefined('x');x = 1:length(y);end
  
% get length of x
n = length(x);

% get y upper and lower bounds
if ~ieNotDefined('yError'),yLow=yError;yHigh = yError;end
if ieNotDefined('yLow'),yLow = zeros(1,n);end
if ieNotDefined('yHigh'),yHigh = yLow;end
if ieNotDefined('yErrorBarType') yErrorBarType = 'both';end

% get x upper and lower bounds
if ~ieNotDefined('xError'),xLow=xError;xHigh = xError;end
if ieNotDefined('xLow'),xLow = zeros(1,n);end
if ieNotDefined('xHigh'),xHigh = xLow;end

% colors and symbols
if ieNotDefined('Symbol'),Symbol = 'o-';end
if ieNotDefined('Color')
  if ~ieNotDefined('MarkerFaceColor')
    Color = MarkerFaceColor;
  else
    Color = 'k';
  end
end
if ieNotDefined('MarkerEdgeColor'),MarkerEdgeColor=Color;end
if ieNotDefined('MarkerFaceColor'),MarkerFaceColor=Color;end
if ieNotDefined('MarkerSize'),MarkerSize=8;end
if ieNotDefined('LineWidth'),LineWidth = 0.5;end

% whether to draw tees or not
if ieNotDefined('tee'),tee = 0;end
if tee
  if ieNotDefined('yTeelen'),yTeelen = mean(diff(x))/10;end
  if ieNotDefined('xTeelen'),xTeelen = mean(diff(y))/10;end
end
hold on
% plot the y error bars
if any(any(yLow ~= 0)) || any(any(yHigh ~= 0))
  for i = 1:length(x)
    switch yErrorBarType
      case {'both','b'}
       plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
     case {'lower','lo','l','bottom','bot'}
       plot([x(i) x(i)],[y(i)-yLow(i) y(i)],'-','Color',Color,'LineWidth',LineWidth);
     case {'higher','upper','up','top','hi'}
       plot([x(i) x(i)],[y(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
      case {'logy'}
       if (y(i)-yLow(i)) > 0
	 plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
       else
	 disp(sprintf('(myerrorbar) Dropping lower errorbar on %i which goes to %f',i,y(i)-yLow(i)));
	 plot([x(i) x(i)],[y(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
       end	 
     otherwise
       plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
    end      
    % draw the tees if necessary
    if tee
      plot([x(i)-yTeelen/2 x(i)+yTeelen/2],[y(i)-yLow(i) y(i)-yLow(i)],'-','Color',Color,'LineWidth',LineWidth);
      plot([x(i)-yTeelen/2 x(i)+yTeelen/2],[y(i)+yHigh(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
    end
  end
end

% plot the x error bars
if any(any(xLow ~= 0)) || any(any(xHigh ~= 0))
  for i = 1:length(x)
    plot([x(i)-xLow(i) x(i)+xHigh(i)],[y(i) y(i)],'-','Color',Color,'LineWidth',LineWidth);
    % draw the tees if necessary
    if tee
      plot([x(i)-xLow(i) x(i)-xLow(i)],[y(i)-xTeelen/2 y(i)+xTeelen/2],'-','Color',Color,'LineWidth',LineWidth);
      plot([x(i)+xHigh(i) x(i)+xHigh(i)],[y(i)-xTeelen/2 y(i)+xTeelen/2],'-','Color',Color,'LineWidth',LineWidth);
    end
  end
end

% plot the symbols
plot(x,y,Symbol,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'Color',Color,'MarkerSize',MarkerSize,'LineWidth',LineWidth);





