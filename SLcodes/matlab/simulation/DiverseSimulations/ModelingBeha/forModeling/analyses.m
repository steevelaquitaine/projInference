

% analyses.m
%
%        $Id: analyses.m 750 2012-09-08 07:12:46Z steeve $
%      usage: analyses
%         by: steeve laquitaine
%       date: 22/10/12
%  copyright: (c) 2012 steeve laquitaine
%    purpose: plot sample dir against est dirs, figure
% inspired from Jazayeri et al, NN.
%
%--------------------------------------------------------------------------
% input
%--------------------------------------------------------------------------
% data: the nm of the file contains the parameters of the data collection
% '120908_Steeve_exp01_data_run1_Pstd1000_mean140_coh024_t250_sess1'

% fig.nm: is a copy of the nm of the file in which 'resul' replace 'data'.
% '120910_Steeve_exp02_data_run8_Pstd1000_mean140_coh024_t250_sess2'


%--------------------------------------------------------------------------
% process
%--------------------------------------------------------------------------
% organize all data into one matrix of cells.
% run specific analyses on the data
% store the results
% publish the project as a .html report

%--------------------------------------------------------------------------
% usage
%--------------------------------------------------------------------------
% % if data are specified
%     data='120910_Steeve_exp02_data_run8_Pstd1000_mean140_coh024_t250_sess2';
%     fig.nm='120910_Steeve_exp02_resul_run8_Pstd1000_mean140_coh024_t250_sess2';
%     analyses(data,fig)

%     % if data are not specified
%     fig.nm='120910_Steeve_resul_databank';
%     analyses(fig)

%--------------------------------------------------------------------------
% issues
%--------------------------------------------------------------------------
% * check what est dir is registered when I choose 0 degree: 0 or
% 360?
% * I don't have access to getang on my laptop even though svn is
% installed ? I wrote an equivalent function getangle
% * think of what plot is good to 1) describe the data and visualize the
% 2) answer to our question.

%--------------------------------------------------------------------------
% abbreviations used in the code
%--------------------------------------------------------------------------
% see http://www.all-acronyms.com/reverse/LVLS
% FA    : factors
% LVLS  : levels
% T     : trial
% est   : estimate

%--------------------------------------------------------------------------
% Stuff to do
%--------------------------------------------------------------------------
% - code to visualize raw data
% - add code lines to check if the FAs are correctly input(orthograph etc..)
% - write function 'graphStat'
% - allow to input only one or two FAs (not always three)

% - structure databank to trash away the dir.series size problem.


%%  Call analyses
function [fig, GenStat] = analyses(data,FAs,fig)
% Check arguments
% if ~any(nargin == [2])
%     help analyses
%     return
% end

global databank
% Gather all available data in directory.
fprintf('% \n %s \n','Now creating a databank...')
Makedatabank(data);

% % Draw estimated dirs = f(displayed dirs)
% fprintf('% \n %s \n','Now identifying the factors"...')
[fig, GenStat]  = initFAs(FAs,fig);

% Describe the patterns in the data.
fprintf('% \n %s \n','Now Describing the patterns in the data...')
DescribeData(databank)

% % Look for "sequential effects"
% fprintf('% \n %s \n','Now testing for "sequential effect"...')
% SequentialEffect







%% Create a matrix of cells in which columns are the variables (databank)
function Makedatabank(data)
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

databank.data =[];

%% Loop over the subfolders to analyse
for j = 1: Folds.nb
    
    % directory of the subfolder to open
    Folds.path(j)= strcat(pathFold,'/',Folds.name(j));
    
    % switch to the subfolder to open
    cd(Folds.path{j})
    
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
    sample_dir_coord = {};
    coh = {};
    
    % check if data have been specified in the function argin
    if ~isempty(datatmp.nm)
        datalisting = datatmp.nm;
        %     datalisting = {datalisting};
        %     display(["--- A databank is being created with the data specified ----']);
        
        % if data have not been specified, gather data from directory
    else
        datalisting = dir('*data*.mat');
        datalisting = {datalisting.name};
        %     display(['--- No data have been specified. A databank is being created with all data in the directory and analysed ---']);
    end
    
    % tic
    
    %%
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
        
        % get the displayed dirs (group 1)
        % when design is balanced across coherences (same trial number).
        if isfield(datai{1,2}.parameter,'dir')
            sample_dir_i = num2cell(datai{1,2}.parameter.dir)';
            sample_dir   = [sample_dir; sample_dir_i];
            % when design is adjusted according to coherences (different trial nb).
        elseif isfield(datai{1,2}.randVars,'myRandomDir')
            sample_dir_i = num2cell(datai{1,2}.randVars.myRandomDir)';
            sample_dir   = [sample_dir; sample_dir_i];
        end
        
        
        % Convert the displayed dirs from degrees to cartesian
        % coordinates
        r = 2.5; % radius of the circular patch
        for jk = 1: size(sample_dir,1)
            sample_dir_coord{jk,1} = num2cell(polar2cartesian(sample_dir{jk},r));
        end
               
        
        
        % Get the displayed coherences (group 2)
        % when design is balanced across coherences (same trial number).
        if isfield(datai{1,2}.parameter,'coherence')
            coh_i = num2cell(datai{1,2}.parameter.coherence)';
            coh  = [coh; coh_i];
            % when design is adjusted according to coherences (different trial nb).
        elseif isfield(datai{1,2}.randVars,'myRandomCoh')
            coh_i = num2cell(datai{1,2}.randVars.myRandomCoh)';
            coh  = [coh; coh_i];
        end
        %%
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
        {'est_dir'      },...
        {'sample_dir_coor'}];
    
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
        est_dirf4Visu   ...
        sample_dir_coord];
    
    % store data over subjects
    databank.data = [databank.data; databanki.data];
end

% Data are organized and saved in a file called datbank in the directory
% switch back to the mother directory
cd(pathFold)

% backup the databank
save ('datbank','databank');

%% Set the factors (predictors) of the analysis
function [fig, GenStat] = initFAs(FAs,fig)
global databank
% Get est data (e.g., estimated dirs)
est_data    = cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));

% Set FAs
% FA 1
FA.g1.thisT = cell2mat(databank.data(:,(strcmp(databank.nm,FAs{1}  ))==1));
FA.g1.nm        = FAs{1};
% FA 2
FA.g2.thisT = cell2mat(databank.data(:,(strcmp(databank.nm,FAs{2}   ))==1));
FA.g2.nm        = FAs{2};
% FA 3
FA.g3.thisT = cell2mat(databank.data(:,(strcmp(databank.nm,FAs{3}  ))==1));
FA.g3.nm        = FAs{3};


% Get and order levels of each group
% group 1 (e.g., priors)
FA.g1.lvlsnm=unique(FA.g1.thisT); %names
FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend'); %order
FA.g1.lvlsnb=numel(FA.g1.lvlsnm);   %number
clear i
for i=1:FA.g1.lvlsnb
    index.g1.lvl_i(i) = {find(FA.g1.thisT==FA.g1.lvlsnm(i))};
end
% group 2 (e.g., coherences)
FA.g2.lvlsnm=unique(FA.g2.thisT); %names
FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend'); %order
FA.g2.lvlsnb=numel(FA.g2.lvlsnm);   %number
for i=1:FA.g2.lvlsnb
    index.g2.lvl_i(i) = {find(FA.g2.thisT==FA.g2.lvlsnm(i))};
end
% group 3 (e.g., displayed dirs)
FA.g3.lvlsnm=unique(FA.g3.thisT); %names
FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'ascend'); %order
FA.g3.lvlsnb=numel(FA.g3.lvlsnm);   %number
for i=1:FA.g3.lvlsnb
    index.g3.lvl_i(i) = {find(FA.g3.thisT==FA.g3.lvlsnm(i))};
end


% Calculate coordinates of average estimated dirs for each condition
% organize data in following order: group 1(subplots) - group 2(colors) -
% group 3(x-axis). Each cell contains repetitions of a condition.
for k=1:FA.g1.lvlsnb
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            index.g1g2g3(j,i,k)={intersect( intersect( index.g1.lvl_i{k},index.g2.lvl_i{j} ),index.g3.lvl_i{i} )};
            
            % calculate statistics of estimated dirs (j:group 2, i:group 3, k:group 1)
            est_dir{j,i,k} = statcircular(est_data(index.g1g2g3{j,i,k},:));
            
        end
    end
end
[fig,GenStat] = graphBasics(FA,est_dir,fig);

%% Plot the factors' effect on the data
function [figs, GenStat] = graphBasics(FA, est_dir,fig)
global databank
% Draw group1-lvl1 (e.g., prior std=80)
% Enlarge figure for good quality publication
% fig1.hdle=figure('Position', [0 0 1000 400]); % pixels
fig1.hdle=figure('color','w'); % pixels
fig1.nm=[fig.nm,'_mean'];

% initialize graphic parameters
% set the colors
c=colormap;
% FA.g2.color = {c(64,:),c(55,:),c(50,:),...
%     c(45,:),c(40,:),c(32,:),...
%     c(27,:),c(22,:),c(15,:),...
%     c(5,:),c(1,:),[0 0 0]}; % group 2(e.g., coherences)
FA.g2.color  = {[0.5 0 0],...
    c(55,:),...
    [1 0.4 0],...
    [1 0.7 0],...
    c(40,:),...
    c(32,:),...
    c(27,:),...
    c(22,:),...
    c(8 ,:),...
    c(5 ,:),...
    c(1 ,:)};

% check if there is enough colors
if numel(FA.g2.color) < FA.g2.lvlsnb
    disp(['--- You may want to add ',...
        num2str(FA.g2.lvlsnb - numel(FA.g2.color)),...
        'more colors to the color code ---']);
    return
end
%% ----------------------------------------------------------------------
% Plot data - mean
%----------------------------------------------------------------------
display('Now plotting the mean of the estimated data...')
% set axes positions
width=1/(FA.g1.lvlsnb+1);
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% Loop over the factors
for k = 1 : FA.g1.lvlsnb %e.g., prior' strength
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    hold all
    % draw priors' mean
    p1 = plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],[databank.data{1,7} databank.data{1,7}],...
        'b--',...
        'linewidth', 1);
    
    % draw cardinal dir
    p2 = plot([0 30   ],[0   0  ],'color',[0.85 0.85 0.85],'linestyle','--');
    p3 = plot([60 120 ],[90  90 ],'color',[0.85 0.85 0.85],'linestyle','--');
    p4 = plot([150 210],[180 180],'color',[0.85 0.85 0.85],'linestyle','--');
    p5 = plot([240 300],[270 270],'color',[0    0    0.85],'linestyle','--',...
        'linewidth', 1);
    
    % add linecode here for drawing an arrow indicating initial position
    % of response line.
    
    % draw ideal performances
    p10 = plot(FA.g3.lvlsnm,FA.g3.lvlsnm','k:');
    
    % organize data
    % ----------------------------------------------------------------
    for j = 1 : FA.g2.lvlsnb %e.g.,coherences
        % loop over  the levels of factor 3
        for i = 1 : FA.g3.lvlsnb %displayed dirs
            % collect data to plot
            fig1.datamean.deg(j,i,k) = est_dir{j,i,k}.deg.mean;
            fig1.datasem (j,i,k) = est_dir{j,i,k}.deg.sem;
            fig1.datamean.coord{j,i,k} = est_dir{j,i,k}.coord.mean;
            
            % extract information about data
            % g1
            fig1.dataInfog1{j,i,k} = FA.g1.lvlsnm(k);
            % g2
            fig1.dataInfog2{j,i,k} = FA.g2.lvlsnm(j);
            % g3
            fig1.dataInfog3{j,i,k} = FA.g3.lvlsnm(i);
            
            % Collect displayed dirs in cartesian coordinates
            r = 2.5;
            fig1.displayed.coord{j,i,k} = polar2cartesian(fig1.dataInfog3{j,i,k},r);
            
            
            % Normalize (linearize) data for plotting (because circular data)
            % ----------------------------------------------------------------\
            % Algorithmic method
            %         % In case the dir is displayed in the 3rd quarter.
            %         if fig1.dataInfog3{j, i, k} > 270 && fig1.datamean.deg(j, i, k)  < 90
            %             % "Linearize" the value of the estimated dir
            %             % fig1.datamean.deg{j, i, k} = 360 + fig1.datamean.deg(j, i, k) ;
            %             fig1.datamean.deg(j, i, k) = 360 + fig1.datamean.deg(j, i, k) ;
            %             % In case the dir is displayed in the 1st quarter.
            %         elseif fig1.dataInfog3{j, i, k} < 90 && fig1.datamean.deg(j, i, k)  > 270
            %             % "Linearize" the value of the estimated dir
            %             %fig1.datamean.deg{j, i, k} = fig1.datamean.deg(j, i, k)  - 360;
            %             fig1.datamean.degtmp(j, i, k) = fig1.datamean.deg(j, i, k)  - 360;
            %         end
            
            % Geometric method
            % Calculate distance between displayed and estimated dirs.
            fig1.Disp_dataDistance{j,i,k} = vectors2signedAngle(fig1.datamean.coord{j,i,k}, fig1.displayed.coord{j,i,k} );
            % Linearize
            fig1.datamean.degLinear(j,i,k) = fig1.dataInfog3{j,i,k} + fig1.Disp_dataDistance{j,i,k};
            
            % Select linearized data or raw data 
%             % linearized
%             fig1.datamean.deg(j,i,k) = fig1.datamean.degLinear(j,i,k);
            %raw data
            fig1.datamean.deg(j,i,k) = fig1.datamean.deg(j,i,k);
        end
        
        % Linear fit
        % ----------------------------------------------------------------
        %%Plot fit for circular data (wrong !)
%         linefit(FA.g3.lvlsnm',fig1.datamean.deg(j,:,k),FA.g2.color{j});
        % Plot fit for linearized data 
        linefit(FA.g3.lvlsnm',fig1.datamean.degLinear(j,:,k),FA.g2.color{j});
        
   
        % Plot data (e.g., mean estimated dir)
        % ----------------------------------------------------------------
        % note: errorbars are sem.
        
        %%Plot data 
%         p12(j) = errorbar(FA.g3.lvlsnm, fig1.datamean.deg(j,:,k), fig1.datasem(j,:,k),'-o',... %groups 3 forms x-axis
%             'color',FA.g2.color{j},...
%             'markerfacecolor',FA.g2.color{j},...
%             'markersize',15  ,...
%             'MarkerEdgeColor','w',...
%             'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        
        myerrorbar(FA.g3.lvlsnm, fig1.datamean.deg(j,:,k), 'yError', fig1.datasem(j,:,k),'Symbol=o-',... %groups 3 forms x-axis
            ['Color=',num2str(FA.g2.color{j})],...
            'LineWidth=0.5',...
            ['MarkerFacecolor=',char(FA.g2.color{j})],...
            'MarkerSize=6',...
            'MarkerEdgeColor=w')%,...
            %'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        
        % legend
        %leg=legend (p12,'location','northwest');
        %leg=legend ('location','northwest');
        %set(leg, 'Box', 'off');
    end
    
    % set graphics
    % set x-axis
    xunit=1:6:FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'ytick',FA.g3.lvlsnm(xunit),'yticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    
    % set labels
    if k==1
        ylabel(ax(k),'Estimated direction (?)','fontsize',14);
    end
    xlabel('Displayed direction (?)','fontsize',14);
    drawPublishAxis
end

% title
ax(k+1,:) = axes('position',[0 0 1 1],'visible','off');
for k = 1 : FA.g1.lvlsnb
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
        'fontsize',20);
    
end

%% ----------------------------------------------------------------------
% Plot data - std
%----------------------------------------------------------------------
display('Now plotting the std of the estimated data...')

% Set figure
fig2.hdle = figure('color','w');
fig2.nm = [fig.nm, '_std'];

% Set the axes positions
width = 1/(FA.g1.lvlsnb+1);

% Arrange the axes
gap = (1-(FA.g1.lvlsnb*width))/(FA.g1.lvlsnb+1);
for k=1:FA.g1.lvlsnb
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

% Set the number of bootstrap
numboot = 100;
tic

% .......
% Draw
% .......
% Loop over the factors
for k = 1 : FA.g1.lvlsnb %e.g., prior' strength
    ax(k) = axes('position',axs.position(k,:));
    axis square
    
    hold all
    % organize the data for plotting
    % loop over the levels of factor 2 (e.g.,coherences)
    for j = 1 : FA.g2.lvlsnb 
        % loop over the levels of factor 3 (e.g.,displayed dirs)
        for i = 1 : FA.g3.lvlsnb 
            % extract information about data
            fig2.datastd (j,i,k) = est_dir{j,i,k}.deg.std;
            fig2.datanb  (j,i,k) = est_dir{j,i,k}.num;                     
            % bootstrap the data to get a distribution of std
            % initialize the data to bootstrap
            data2boot = [];
            fig2.databootstdSeries{j,i,k} = [];
            % collect the data to bootstrap (e.g. data at each displayed dir)
            fig2.datacoord (j,i,k) = est_dir{j,i,k}.coord;
            data2boot.data = fig2.datacoord (j,i,k).all;
            % create samples
            for jj = 1 : numboot
                % count the data
                data2boot.nb = size(data2boot.data,1);
                % if data are observed in x and sample >= 2
                if isempty(data2boot.data) == 0 && data2boot.nb >= 2     
                    % bootstrap the data in x
                    % pick randomly one data from x 'data2boot.nb' times
                    for iii = 1 : data2boot.nb
                        jit = ceil(rand * data2boot.nb);
                        data_boottmp{j,i,k}{jj}(iii,:) = data2boot.data(jit,:);
                        % store for fig2
                        fig2.databoot{j,i,k}{jj}(iii,:) = data_boottmp{j,i,k}{jj}(iii,:) ;
                        fig2.databootnb{j,i,k}{jj}(iii,:) = data2boot.nb;
                    end
                    % calculate the std for each bootstrap
                    fig2.databootstd{j,i,k}{jj} = statcircular(fig2.databoot{j,i,k}{jj});
                    % collect the std in a vector
                    fig2.databootstdSeries{j,i,k} = [fig2.databootstdSeries{j,i,k} fig2.databootstd{j,i,k}{jj}.deg.std];
                    % if there are no data in x
                else
                    % fill in a nan
                    data_boottmp{j,i,k}{jj} = nan;
                    % store for fig2
                    fig2.databoot{j,i,k}{jj} = data_boottmp{j,i,k}{jj} ;
                    % there are no data
                    fig2.databootnb{j,i,k}{jj} = 0 ;
                    % nan the std
                    fig2.databootstd{j,i,k}{jj} = nan;
                    % nan the vector of std
                    fig2.databootstdSeries{j,i,k} = nan;
                end
            end
            % calculate the std over bootstrapped samples
            fig2.databootMeanOfstd{j,i,k} = nanmean(fig2.databootstdSeries{j,i,k});
            fig2.databootStdOfstd{j,i,k} = std(fig2.databootstdSeries{j,i,k});
            
            %% ---- testing
            
            % extract information about data
            % g1
            fig2.dataInfog1{j,i,k} = FA.g1.lvlsnm(k);
            % g2
            fig2.dataInfog2{j,i,k} = FA.g2.lvlsnm(j);
            % g3
            fig2.dataInfog3{j,i,k} = FA.g3.lvlsnm(i);
        end
        % plot data (e.g., std)
        % bootstrap to get errorbars.
        
%         p12(j) = errorbar(FA.g3.lvlsnm,[fig2.databootMeanOfstd{j,:,k}], [fig2.databootStdOfstd{j,:,k}],...
%             '-o',... %groups 3 forms x-axis
%             'color', FA.g2.color{j},...
%             'markerfacecolor', FA.g2.color{j},...
%             'markersize', 15,...
%             'MarkerEdgeColor','w',...
%             'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));
        
        myerrorbar(FA.g3.lvlsnm,[fig2.databootMeanOfstd{j,:,k}], 'yError',[fig2.databootStdOfstd{j,:,k}],...
            'Symbol=o-',...
            'LineWidth=0.5',...
            ['Color=',num2str(FA.g2.color{j})],...
            ['MarkerFacecolor=', num2str(FA.g2.color{j})],...
            'MarkerSize=6',...
            'MarkerEdgeColor=w')%,...
            %'displayname', strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))));

        %leg = legend (p12,'location','northwest');
        %leg = legend ('location','northwest');
        %set(leg, 'Box', 'off');
    end
    
    %set graphics
    %set x-axis
    xunit = 1: 6: FA.g3.lvlsnb;
    set(gca,...
        'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
        'fontsize',20);
    %limits
    xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);
    ylim([-1 200]);
    
    % set labels of axes
    if k==1
        ylabel(ax(k),'std of estimated direction','fontsize',14);
    end
    xlabel('Displayed direction (?)','fontsize',14);
    drawPublishAxis
end
% toc

% title
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
        'fontsize',20);
end

% Draw raw data
display('Now plotting the mean estimated data...')
[fig4] = drawRawEstDir(fig2);

% backup figures informations
figs = {fig1, fig2};
% %%
% %----------------------------------------------------------------------
% % Stat of ANCOVA: effect of FA1(e.g., coherence) on the slope of data=f(FA3)).
% %----------------------------------------------------------------------
% % Set x: (e.g.,FA3 - displayed dirs)
% x = [];
% x = reshape(repmat(FA.g3.lvlsnm , FA.g2.lvlsnb*FA.g1.lvlsnb,1), [],1);
% 
% % Set y: (data; e.g., estimated dirs)
% % Loop over the levels of Factor 1
% y = [];
% for k = 1 : FA.g1.lvlsnb 
%     ytmp = [];
%     ytmp = reshape(fig1.datamean.deg(:,:,k)',[],1);
%     y = [y; ytmp];
% end
% 
% % Set g the groups (e.g., FA1 - coherences)
% g = [];
% g = reshape(repmat(1: FA.g1.lvlsnb, FA.g3.lvlsnb*FA.g2.lvlsnb,1), [],1);
% 
% % Run the Ancova
% Ancova_effect_of_FA1 = statAncova(x, y, g);
% display('Your Ancova is being performed...')
% 
% % % Report the statistics of the Ancova
% % % Degrees of freedom
% % % Within groups  (e.g., = f(nb of subjects))
% % dfwgT1   = GenStat{2,1}{2,2}{5,2};
% % % Between groups  (e.g.,= f(nb of group))
% % dfbgT1   = GenStat{2,1}{2,2}{4,2};
% % % F-value for the interaction effect
% % FvalueT1 = GenStat{2,1}{2,2}{4,5};
% % % p-value of the tested interaction effect
% % pvalueT1 = GenStat{2,1}{2,2}{4,6};
% % GenStat2ReportT1 = {dfbgT1, dfwgT1, FvalueT1, pvalueT1};
% % display(GenStat2ReportT1)

%%
%----------------------------------------------------------------------
% Stat of ANCOVA: effect of FA1(e.g., prior's width) on the slope of data=f(FA3)).
%----------------------------------------------------------------------
% Set x: (e.g.,FA3 - displayed dirs)
x = [];
x = reshape(repmat(FA.g3.lvlsnm , FA.g2.lvlsnb*FA.g1.lvlsnb,1), [],1);

% Set y: (data; e.g., estimated dirs)
% Loop over the levels of Factor 1
y = [];
for k = 1 : FA.g1.lvlsnb 
    ytmp = [];
    ytmp = reshape(fig1.datamean.deg(:,:,k)',[],1);
    y = [y; ytmp];
end

% Set g the groups (e.g., FA1 - coherences)
g = [];
g = reshape(repmat(1: FA.g2.lvlsnb, FA.g3.lvlsnb*FA.g1.lvlsnb,1), [],1);

% Run the Ancova
Ancova_effect_of_FA2 = statAncova(x, y, g);
display('Your Ancova is being performed...')

% % Report the statistics of the Ancova
% % Degrees of freedom
% % Within groups  (e.g., = f(nb of subjects))
% dfwgT1   = GenStat{2,1}{2,2}{5,2};
% % Between groups  (e.g.,= f(nb of group))
% dfbgT1   = GenStat{2,1}{2,2}{4,2};
% % F-value for the interaction effect
% FvalueT1 = GenStat{2,1}{2,2}{4,5};
% % p-value of the tested interaction effect
% pvalueT1 = GenStat{2,1}{2,2}{4,6};
% GenStat2ReportT1 = {dfbgT1, dfwgT1, FvalueT1, pvalueT1};
% display(GenStat2ReportT1)


%% ----------------------------------------------------------------------
% Stat T1: ANCOVA (effect of FA2 on the slope of data=f(FA3) compared to 1.
%----------------------------------------------------------------------
% e.g., FA 1 is coherence
% e.g., FA 2 is prior's strength
% e.g., FA 3 is displayed dir

% Loop over each condition FA1-FA2
% Levels of FA 1
for k = 1 : FA.g1.lvlsnb
    % Levels of FA 2
    for j = 1 : FA.g2.lvlsnb
    % Set x: displayed dirs (levels of FA3)
    x = [];
    x = reshape(repmat(FA.g3.lvlsnm, 2, 1), [],1);
    % Set y: data (e,g., estimated dirs and dirs with slope 1).
    y = [];
    y = reshape([fig1.datamean.deg(j,:,k)', FA.g3.lvlsnm], [],1);
    % Set the groups (g)
    g = [];
    g = reshape(repmat([1,2], FA.g3.lvlsnb,1), [],1);
%     % Title rows
%     GenStatAgainst1{j+1,  1} = [FA.g2.nm num2str(FA.g2.lvlsnm(j))];
%     % title columns
%     GenStatAgainst1{  1,k+1} = [FA.g1.nm num2str(FA.g1.lvlsnm(k))];
    % run analysis
    Ancova_effect_of_FA2vs1{j+1,k+1} = statAncova(x,y,g);
    end
end

% Stat to report (write sthg here)


%% ---------------------------------------------------------------------
% Stat: ANCOVA - compare the slopes of FA 2's levels for each FA 1.
%---------------------------------------------------------------------
% e.g., FA 1 is coherence
% e.g., FA 2 is prior's strength
% e.g., FA 3 is displayed dir
% Loop over the levels of FA1 (e.g., coherences)
for k = 1 : FA.g1.lvlsnb
    % Set x: FA3 - e.g., displayed dirs (levels of FA 3)
    x = [];
    x = reshape(repmat(FA.g3.lvlsnm, FA.g2.lvlsnb,1),[],1);
    
    % Set y: data e.g., estimated dirs
    y = [];
    y = reshape(fig1.datamean.deg(:,:,k)', [], 1);
    
    % Set g: groups (levels of FA2, e.g., priors)
    g = [];
    g = reshape(repmat(FA.g2.lvlsnm', FA.g3.lvlsnb,1), [],1);
    
    % title columns
    GenStatBetwLv{1,k} = [FA.g1.nm num2str(FA.g1.lvlsnm(k))];
    % run analysis
    GenStatBetwLv{2,k} = statAncova(x,y,g);
end

% Stat to report(write sthg here)


%----------------------------------------------------------------------
% Store the statistics
%----------------------------------------------------------------------
% Title columns
GenStat = [ {'Ancova_effect_of_FA2'},...
    {'Ancova_effect_of_FA2vs1'},...
    {'Ancova_lvl1_vs_lvl2_FA1'}];
% Data
GenStat(2,:) = [ {Ancova_effect_of_FA2},...
    {Ancova_effect_of_FA2vs1},...
    {GenStatBetwLv}];

display('The results of the Ancova have been stored')


% Report the statistics of the Ancova
% Degrees of freedom
% Within groups  (e.g., = f(nb of subjects))
dfwgT1   = GenStat{2,1}{2,2}{5,2};
% Between g roups  (e.g.,= f(nb of group))
dfbgT1   = GenStat{2,1}{2,2}{4,2};
% F-value for the interaction effect
FvalueT1 = GenStat{2,1}{2,2}{4,5};
% p-value of the tested interaction effect
pvalueT1 = GenStat{2,1}{2,2}{4,6};
GenStatT1 = {dfbgT1, dfwgT1, FvalueT1, pvalueT1};
display(GenStatT1)


% Report the statistics of the Ancova
% Degrees of freedom
% Within groups  (e.g., = f(nb of subjects))
dfwgT2   = GenStat{2,2}{2,2}{2,2}{5,2};
% Between groups  (e.g.,= f(nb of group))
dfbgT2   = GenStat{2,2}{2,2}{2,2}{4,2};
% F-value for the interaction effect
FvalueT2 = GenStat{2,2}{2,2}{2,2}{4,5};
% p-value of the tested interaction effect
pvalueT2 = GenStat{2,2}{2,2}{2,2}{4,6};
GenStatT2 = {dfbgT2, dfwgT2, FvalueT2, pvalueT2};
display(GenStatT2)


% Report the statistics of the Ancova
% Degrees of freedom
% Within groups  (e.g., = f(nb of subjects))
dfwgT3   = GenStat{2,3}{2,1}{2,2}{5,2};
% Between groups  (e.g.,= f(nb of group))
dfbgT3   = GenStat{2,3}{2,1}{2,2}{4,2};
% F-value for the interaction effect
FvalueT3 = GenStat{2,3}{2,1}{2,2}{4,5};
% p-value of the tested interaction effect
pvalueT3 = GenStat{2,3}{2,1}{2,2}{4,6};
GenStatT3 = {dfbgT3, dfwgT3, FvalueT3, pvalueT3};
display(GenStatT3)

% Under construction here !

%% ----------------------------------------------------------------------
% Backup figure
%----------------------------------------------------------------------
autobackup(fig1.hdle,fig1.nm,'.fig');

% [fig]=graphStat(FA, est_dir,fig);
% Graph the quantitative effect of the prior's strength
if FA.g1.lvlsnb == 2
    [fig3] = graphPriorEffect(fig1,FA);
end

%% Plot the effect of Factor f on the data (e.g., prior's strength)
function [fig3] = graphPriorEffect(fig1,FA)
% Calculate the effect of the first FA on the data (level 2 - level 1)
% note: it is possible two compare the data for two levels only.

% set figure
% fig3.hdle=figure('Position', [0 0 1000 400]);%pixels
fig3.hdle=figure('color','w'); % pixels
axis square
fig3.nm=[fig1.nm,'_PriorEffect'];

% Check if the first FA has exactly 2 levels.
if FA.g1.lvlsnb ~= 2
    disp([ '----- graphPriorEffect has been cancelled. This analysis is valid only if the number of levels in the first FA = 2 ----'])
    return
end
fig3.datamean = [fig1.datamean.deg(:,:,1) - fig1.datamean.deg(:,:,2)]'; %e.g., prior's width effect


% Draw linear fit to the data
hold all;
for j = 1:FA.g2.lvlsnb
    linefit(FA.g3.lvlsnm,fig3.datamean(:,j),FA.g2.color{j});
end

% Draw the prediction for an absence of effect of the first FA (e.g., prior width)
plot(FA.g3.lvlsnm,zeros(FA.g3.lvlsnb,1),...
    '--','color',[0.5 0.5 0.5]);

% Draw the mean of the prior
plot([140 140],[min(min(fig3.datamean)) max(max(fig3.datamean))],'b--',...
    'linewidth', 1);

% Draw the effect of the first FA(e.g., prior's width)
for j=1:FA.g2.lvlsnb
    p13(j)=plot(FA.g3.lvlsnm,fig3.datamean(:,j),...
        'o',...
        'markersize',10  ,... % 18 for illustrator (9 for publishing)
        'markerfacecolor', FA.g2.color{j},...
        'MarkerEdgeColor','w',...
        'color', FA.g2.color{j},...
        'displayname',strcat(FA.g2.nm,'--',num2str(FA.g2.lvlsnm(j))));
end

% Set the legend
leg=legend (p13,'location','northwest');
set(leg, 'Box', 'off');

% Set the axes label
ylabel('Difference of estimated directions (?, narrow - wide)','fontsize',14);
xlabel('Displayed directions (?)','fontsize',14);

% set the graph parameters
% set the unit of the x-axis
xunit=1:6:FA.g3.lvlsnb;
set(gca,...
    'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
    'fontsize',20);

xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);

% backup figures and parameters
autobackup(fig3.hdle,fig3.nm,'.fig');

%% Test the hypothesis of "sequential effect"
function SequentialEffect
global databank

% Get displayed dirs
dir.dispDegr    = cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));

% Set the factors
FAs{1} = 'Pstd';
FAs{2} = 'coh';
% FAs{1} = 'coh';
% FAs{2} = 'Pstd';

% Set FAs
% FA 1
FA.g1.thisT = cell2mat(databank.data(:,(strcmp(databank.nm,FAs{1}  ))==1));
FA.g1.nm        = FAs{1};
FA.g1.unique = unique(FA.g1.thisT);
% FA 2
FA.g2.thisT = cell2mat(databank.data(:,(strcmp(databank.nm,FAs{2}   ))==1));
FA.g2.nm        = FAs{2};
FA.g2.unique = unique(FA.g2.thisT);

% set the radius of the circle
r = 4; % but the true size is 2.5 degrees (just for visualization purpose)
% Get the coordinates of the displayed dirs
dir.dispCoor = polar2cartesian(dir.dispDegr,r);
% Get the coordinates of the most likely dir (prior's mean)
dir.priormeanCoor = polar2cartesian(225,r);
% Get estimated dirs
dir.estiCoor = cell2mat(databank.data(:,5));

% Get dirs displayed between the most likely (prior's mean) and
% the previous dirs displayed further away from the most likely.
% e.g., 225 degrees (most likely)< 180 deg (new) <135 deg (old)
% Initialize a matrix to store the trials of each conditions
ibkp.leftward = []; 
% Loop over the displayed dirs
for i = 2 : size(dir.dispCoor,1)
    %if the new(t) and the old(t-1) dirs are displayed in the
    %hemicircle that is leftward to the most likely dir (prior's mean).
    % i.e., 0< angle(prior, new) <180
    angle.PriortoNew(i) = vectors2signedAngle(dir.priormeanCoor, dir.dispCoor(i,:));
    angle.NewtoOld(i)   = vectors2signedAngle(dir.dispCoor(i,:), dir.dispCoor(i-1,:));
    if  0< angle.PriortoNew(i) && angle.PriortoNew(i) <180 && 0< angle.NewtoOld(i) && angle.NewtoOld(i) <180
        %Store the new trial
        ibkp.leftward = [ibkp.leftward; i];
    end
end


% Look at the sequential effect in different conditions
% Set the factor (e.g., coherence)
FA.sample = FA.g2.unique;
FA.series = FA.g2.thisT;

% Prepare the factor
ibkp.leftward.FAtmp = [];
ibkp.leftward.FA = {};
figure('color','w')

% Loop over the factor
for j = 1 : numel(FA.sample)
    ibkp.leftward.FAtmp = [];
    % Loop over the displayed dirs
    for i = 2 : size(dir.dispCoor,1)
        %if the new(t) and the old(t-1) dirs are displayed in the
        %hemicircle that is leftward to the most likely dir (prior's mean).
        % i.e., 0< angle(prior, new) <180
        angle.PriortoNew(i) = vectors2signedAngle(dir.priormeanCoor,...
            dir.dispCoor(i,:));
        angle.NewtoOld(i)   = vectors2signedAngle(dir.dispCoor(i,:),...
            dir.dispCoor(i-1,:));
        %Select trials when new & old dirs are in the
        %hemicircle leftward to the most likely dir (prior's mean).
        if  0< angle.PriortoNew(i) && angle.PriortoNew(i) <180 && ...
                0< angle.NewtoOld(i) && angle.NewtoOld(i) <180 && ...
                isequal(FA.series(i), FA.sample(j))
            %Store the new trial
            ibkp.leftward.FAtmp = [ibkp.leftward.FAtmp; i];
        end
        
%         %%% TEST the code ---- as predicted so OK.
%         %Select trials when old dirs is the most likely dir (prior's mean).
%         if  angle.PriortoNew(i) == - angle.NewtoOld(i) && ...
%                 isequal(FA.series(i), FA.sample(j))
%             %Store the new trial
%             ibkp.leftward.FAtmp = [ibkp.leftward.FAtmp; i];
%         end
        
        
    end
    %Store the trials for each level of the factor
    ibkp.leftward.FA{j} = ibkp.leftward.FAtmp;
    
    %Calculate the coordinates of the average displayed dir (new),
    %& estimated dir.
    dir.dispNewMean{j} = vectorAverage(dir.dispCoor(ibkp.leftward.FA{j}  ,:));
    dir.dispOldMean{j} = vectorAverage(dir.dispCoor(ibkp.leftward.FA{j}-1,:));
    dir.estiMean{j}    = vectorAverage(dir.estiCoor(ibkp.leftward.FA{j}  ,:));
    
    % Draw the dirs.
    % Intitialize polar
    ax1(j) = subplot(1,numel(FA.sample),j) 
    p(j)= polar(0,r);
    set(p(j), 'Visible', 'Off');
    hold on;
    % Remove unnecessary lines and text.
    % Find the lines in the polar plot
    h = findall(gca,'type','line');
    % Remove the handle for the polar plot line from the array
    h(h == p(j)) = [];
    %     Delete other lines
    delete(h([2,3,5:end]));
    % Find all text objects in the polar plot
    t = findall(gca,'type','text');
    % set the text fontsize
    for i = 1 : numel(t)
        set(t(i),'fontsize',18)
    end
    % Delete text objects
    delete(t(end-1:end));
    
    % set arrows
    arrowsz = 10;
    arrowidth = 3;
    % Predictions are:
    % #Sequential effect (total) - bias toward previous dir.
    % #Bayesian inference - bias toward most likely dir.
    
    % Draw the average displayed dir at t
    % Scale the length of the arrows
    arrowidth = 0;
    arrowsz = 7;
    arrow_length = sqrt(dir.dispNewMean{j}.coord.mean(1)^2 + dir.dispNewMean{j}.coord.mean(2)^2);
    scale2radius = r/arrow_length;
    a(2) = arrow([0,0], scale2radius*[dir.dispNewMean{j}.coord.mean(1), dir.dispNewMean{j}.coord.mean(2)],...
        arrowsz,[],[],...
        arrowidth,...
        'EdgeColor','none',...
        'facecolor',[0.7 0.7 0.7],...
        'displayName','Current dir')
    % Draw the average displayed dir at t-1
    % Scale the length of the arrows
    arrow_length = sqrt(dir.dispOldMean{j}.coord.mean(1)^2 + dir.dispOldMean{j}.coord.mean(2)^2);
    scale2radius = r/arrow_length;
    a(3) = arrow([0,0], scale2radius*[dir.dispOldMean{j}.coord.mean(1), dir.dispOldMean{j}.coord.mean(2)],...
        arrowsz,[],[],...
        arrowidth,...
        'EdgeColor','none',...
        'facecolor',[1 0.7 0.7],...
        'displayName','Previous dir')
    % Draw the most likely dir (prior's mean).
    a(4) = arrow([0,0], [dir.priormeanCoor(1), dir.priormeanCoor(2)],...
        arrowsz,[],[],...
        arrowidth,...
        'EdgeColor','none',...
        'facecolor',[0.7 0.7 1],....
        'displayName','Mean of the prior')
    % Draw the average estimated dir.
    arrow_length = sqrt(dir.estiMean{j}.coord.mean(1)^2 + dir.estiMean{j}.coord.mean(2)^2);
    scale2radius = r/arrow_length;
    a(1) = arrow([0,0], scale2radius*[dir.estiMean{j}.coord.mean(1), dir.estiMean{j}.coord.mean(2)],...
        arrowsz,[],[],...
        arrowidth,...
        'EdgeColor', [0 0 0],...
        'facecolor', [0 0 0],...
        'displayname', 'Estimated dir')
        
    hold off
    %title(num2str(FA.sample(j)),'fontsize',32, 'fontweight','bold')
    text(0,350,num2str(FA.sample(j)))
end
% Legend the graph
legend(a(1:4),'location','BestOutside');
set(gca, 'fontsize',18)
legend boxoff

%% Describe data
function DescribeData(databank)

%% get data(e.g., estimated directions)
dir.es    =cell2mat(databank.data(:,strcmp(databank.nm,'est_dir')==1));
prior.m   =unique(cell2mat(databank.data(:,strcmp(databank.nm,'priormean')==1)));
dir.di  =cell2mat(databank.data(:,strcmp(databank.nm,'sample_dir')==1));
coh       =cell2mat(databank.data(:,strcmp(databank.nm,'coh')==1));
pstd       =cell2mat(databank.data(:,strcmp(databank.nm,'Pstd')==1));

x=5:10:360;

%set colors
cmp=colormap('jet');
colors{1}=cmp;
%color of equidistant direction
colors{2}=cmp(end:-1:1,:);

%Organize data for a clear plot
%The data for diplayed direction that are equidistant to the prior are 
%displayed on the same subplot.
%sort data
C.cohPstd=coh==0.06 & pstd==80;

%range of displayed directions
dir.dicoh=dir.di(C.cohPstd);
dir.di=unique(dir.dicoh);
dir.dicount=numel(dir.di);

%get distances to prior
dist2Pr.sign=dir.di-225;
dist2Pr.abs=abs(dist2Pr.sign);
dist2Pr.uniq=unique(dist2Pr.abs);
dist2Pr.uniqCount=numel(dist2Pr.uniq);

% Initialize figure
scrsz = get(0,'ScreenSize');
figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)])
hold all
for i=1:dist2Pr.uniqCount
    subp(i)=subplot(dist2Pr.uniqCount,1,i);
end
dir.dist2Pr=dir.dicoh-225;

%estimated directions
dir.escoh=dir.es(C.cohPstd);

%Look at each displayed direction
for i=1:dist2Pr.uniqCount
    
    %draw 
    subplot(subp(i))
    
    %directions equidistant to prior are displayed together
    eqpos=find(dist2Pr.abs==dist2Pr.uniq(i));
    for j=1:numel(eqpos)
        
        %sort data
        dir.dispos=dir.dist2Pr==dist2Pr.sign(eqpos(j));
        C.disp=dir.dispos;
        
        %if this condition exists
        if sum(C.disp)~=0;
            
            %calculate pdf....................
            [p xpdf]=makePdf(x,dir.escoh(C.disp),'raw');
            
            %show sample size
            dir.esCount=numel(dir.escoh(C.disp));
            fprintf('\n %12s %12g \n', 'number of trial is', dir.esCount)
           
            %draw.............................
            hold all
            h1=area(xpdf,p,'facecolor',colors{j}(i,:),'edgecolor','none');
            
            %markers
            h2=plot([dir.dicoh(C.disp) dir.dicoh(C.disp)],[0 max(p)],'-','color','k','linewidth',1.5);
            h3=plot([prior.m prior.m],[0 max(p)],'--','color','b','linewidth',1.5);
            
            alpha(0.4)
        end
    end   
    
    %put on top
    uistack(h1, 'top')
    uistack(h2, 'top')
    
    %label    
    ylabel(subp(i),'Probability','fontsize',8)
    if i==dist2Pr.uniqCount
        xlabel('Estimated directions (?)','fontsize',14)
    end
    %drawPublishAxis
end












%% Nested functions
% Backup
function autobackup(task,filename,filetype)
% input:
% structure called 'fig' with 2 fields:
% - 'name': e.g., '120910_Steeve_exp02_resul_run8_var1000_mean195_coh024_t250_sess2'
% - 'hdle': the handle of the figure you want to backup

% check if the file name exists
r=dir;
clear i
nametocheck=strcat(filename,filetype);
for i=1:length(r); scanres(i) = strcmp(r(i).name,nametocheck); end

% if the file name exists, increment
i=0;
if ~isempty(find(scanres==1))
    while ~isempty(find(scanres==1))
        i=i+1;
        % increment the name
        filename=strcat(filename,'_0',num2str(i));
        nametocheck=strcat(filename,filetype);
        % check if the name exists already
        for j=1:length(r); scanres(j) = strcmp(r(j).name,nametocheck); end
    end
    errorms=[' "This filename exists already. The new filename is "',filename,'" '];
    disp(errorms)
end
if strcmp(filetype,'.fig')
    saveas(task,filename,'fig'); %name is a string
elseif strcmp(filetype,'.mat')
    save(f0ilename,inputname(1)); %name is a string
end

% Calculate the circular statistics of the data
function [data] = statcircular(coord)
% input a vector of cartesian coordinates (coord)
% output :
%   -

% register the coordinates of the input dirs
data.coord.all = coord;

% convert from cartesian coordinates to angles (in degree)
data.deg.all = getangle(coord(:,1), coord(:,2));

% calculate the cartesian coordinates of the mean estimated dir.
data.coord.mean(:,1) = nanmean(coord(:,1)); % cartesian coordinate of the mean dir est
data.coord.mean(:,2) = nanmean(coord(:,2));

% calculate the mean dir estimate (in degree)
data.deg.mean = getangle(data.coord.mean(:,1),data.coord.mean(:,2)); %(in degree)


% calculate the std to the mean dir est (in degree); !!! could be a
% subfunction itself
% Apply the rule of thumb that follows. It seems to work fine intuitively. It would be nice to
% fine a cleaner way to calculate the std.
% initialize the 'sample' and 'mean' variables used to calculate the std
% count the number of data
data.num = numel(data.deg.all); % sample size
% collect data for the calculation of the std
data.deg.allforstd = data.deg.all;
% collect the mean of the data for the calculation of the std
data.deg.meanforstd = repmat(data.deg.mean,data.num,1);

% if the resulting mean dir is between 0 and 180.
if data.deg.mean + 180 <= 360
    % sample each est dir
    for i=1:data.num
        % if the est dir sampled is >= mean dir + 180
        if data.deg.all(i) >= data.deg.mean + 180
            data.deg.allforstd(i)=data.deg.all(i) - 360;
        end
    end
    % if the resulting mean dir is between 180 and 360.
else
    % sample each est dir
    for i=1:data.num
        % if the est dir sampled is <= the mean dir - 180
        if data.deg.all(i) <= data.deg.mean-180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

% now calculate the variance of the estimated dir.
data.deg.var = nanmean((data.deg.allforstd - data.deg.meanforstd).^2,1);

% and now calculate the std
data.deg.std = sqrt(data.deg.var);

% and now calculate the sem
data.deg.sem = data.deg.std/sqrt(data.num);

% Convert from cartesian coordinates to angles (degree)
function [angle] = getangle(x,y)
% check! to check if the function works fine, write
% e.g., input=180; output=getangle(cos(angle*pi/180),sin(angle*pi/180));
% if the function works input and output should always be the same between
% 0 and 360.

% convert from cartesian coordinates to angle in radians
angle = atan(y./x); % theta = arctan(opposite/adjacent);

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

% Evaluate the distance between two vectors
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

% Calculate an vector average
function [data] = vectorAverage(coord)
% Calculate the cartesian coordinates of the mean estimated dir.
data.coord.mean = nanmean(coord,1); % cartesian coordinate of the mean dir est

% Calculate the mean dir estimate (in degree)
x = data.coord.mean(:,1);
y = data.coord.mean(:,2);

% radian
data.radian.mean = atan(y./x); %theta = arctan(opposite/adjacent);
% degrees
data.degree.mean = data.radian.mean.* 180/pi;

% Convert from radian to degrees
function [degree] = radian2degree(radians)
degree = radian*pi/180;

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
    'linewidth',1,...
    'color',linecolor);

% Do analysis of covariance (ancova)
function [Stat] = statAncova(x,y,group)
%
% % Input:
% % x: e.g.,displayed dir
% % y: e.g.,estimated dir.
% % group: "strong prior" vs "weak prior".
%
% % OUTPUT: interaction p-value in Anova table specifies if slopes differ.
%
% % Important issue! I can't find how to apply ancova on circular data.
% % Solution: Ancova can be performed with the vector averaged displayed estimated
% % dirs and the displayed dirs.
%
% % x
% databank.data(:,8)
%
% % y
%
% % group
% databank.data(:,6)
%
[h,atab,ctab,stats] = aoctool(x,y,group,0.05,'','','','off');
Stat = {'h','anovatab','coeftab','stats'};
Stat(2,:) = {h,atab,ctab,stats};

% Draw the raw data (e.g., individual dir estimates)
function [fig4] = drawRawEstDir(fig2)

% still under developpment

% input
% angles = vector of cartesian coordinates on the unit circle.
% r = radius of the motion patch: 2.5 degrees.

% Title the panels
figTitles = {'coh32', 'coh16','coh8'};

% set the factors
% factor 1(e.g., coherence)
FA1.lvlnum = size(fig2.datacoord, 1);
% factor 2(e.g., displayed dir)
FA2.lvlnum = size(fig2.datacoord, 2);
% factor 3(e.g., prior)
FA3.lvlnum = size(fig2.datacoord, 3); 

% Initialize the figure
% loop over factor 2 (e.g., coherence)
for  k = 1 : FA3.lvlnum
    fig4(k) = figure('color','w');
    % loop over factor 3 (e.g., displayed dirs)
    for i = 1 : FA2.lvlnum 
        r = 2.5;
        %initialize figure
        subplot(6,6,i)
 
        %intitialize polar
        p = polar(0,r);
        set(p, 'Visible', 'Off');
        
        % %     Remove unnecessary lines and text.
        %     % Find the lines in the polar plot
        %     h = findall(gcf,'type','line');
        %     % Remove the handle for the polar plot line from the array
        %     h(h == p) = [];
        % %     Delete other lines
        %     delete(h);
        %     % Find all text objects in the polar plot
        %     t = findall(gcf,'type','text');
        %     % Delete text objects
        %     delete(t);
        
        % set arrows
        arrowsz = 2;
        arrowidth = 1;
        
        %         FA1.lvlnum = size(figOut{2}.datacoord, 1);
        %         FA2.lvlnum = size(figOut{2}.datacoord, 2);
        %         FA3.lvlnum = size(figOut{2}.datacoord, 3);
        
        % set the colors of factor 3
        c = colormap;
        FA3.lvlcolor = {[0.5 0 0],...
            c(55,:),...
            [1 0.4 0],...
            [1 0.7 0],...
            c(40,:),...
            c(32,:),...
            c(27,:),...
            c(22,:),...
            c(8 ,:),...
            c(5 ,:),...
            c(1 ,:)};
        
        % set parameters for the graphics
        a = [1.3; 1.6; 1.9; 2.1; 2.4; 2.7];
        
        % loop over the levels of factor 3 (e.g., coherence)
        for j = 1 : FA1.lvlnum
            cartCoord.all  = fig2.datacoord(j,i,k).all;
            cartCoord.mean = fig2.datacoord(j,i,k).mean;
            
            arrow_lengthOld = sqrt(cartCoord.mean(1)^2 + cartCoord.mean(2)^2);
            scale2radius = 2.5/arrow_lengthOld;
            
            % draw raw and mean estimated dirs
            if ~isempty(cartCoord.all)
                arrow(0.9*cartCoord.all*a(j),cartCoord.all*a(j), arrowsz, [], [], arrowidth,'EdgeColor', FA3.lvlcolor{j},'facecolor', FA3.lvlcolor{j})
                arrow([0,0],scale2radius*cartCoord.mean, 2, [], [], 1,'EdgeColor', FA3.lvlcolor{j},'facecolor', FA3.lvlcolor{j})
            end
            % draw displayed dirs
            hold on;
            h2 = polar(fig2.dataInfog3{j,i,k}*pi/180, 2*r, 'gs');
            set(h2,'markerfacecolor',[0 0.65 0], 'Markersize',10)
%             leg=legend (p12,'location','northwest');
%             set(leg, 'Box', 'off');
        end
    end
    title (figTitles{k})
    
end    

% Calculate the smoothed pdf of a dataset
function [p xpdf]=makePdf(x,dataset,type)

if strcmp(type,'raw')
    % create pdf (raw)
    c=histc(dataset,x);
    p=c/sum(c);
    xpdf=x;
end

if strcmp(type,'smooth')
    % create pdf (smoothed)
    xpdf=min(x):1:max(x);
    p=ksdensity(dataset,xpdf);
    p=p/sum(p);
end



