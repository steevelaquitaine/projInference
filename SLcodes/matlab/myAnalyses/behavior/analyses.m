%analyses.m
%
%       $Id: analyses.m 750 2012-09-08 07:12:46Z steeve $
%     usage: analyses
%        by: steeve laquitaine
%      date: 22/10/12 last update: 150715
% copyright: (c) 2012 steeve laquitaine
%
%   purpose: run a set of analyses on motion direction
%            estimation data.
%            Data parameters and directories should have been initialized
%            with "SLinitAnalyses" code.
%
%
%
%     usage:
%           ex.
%               %load stim file and run
%               analyses
%
%           ex.
%               SLinitAnalyses
%           ex.
%               analyses({'sub01','sub02'},{'StimStrength','Pstd','FeatureSample'},'experiment','vonMisesPrior','dataDis','signDistance');
%           ex.
%               analyses({'sub02'},{'StimStrength','Pstd','FeatureSample'},'experiment','vonMisesPrior','dataDis','signDistance','dataMeanAndStd');
%           ex.
%               analyses({'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},{'StimStrength','Pstd','FeatureSample'},'experiment','vonMisesPrior','dataDis','absDistance','fastrun');
%   
%           
%          
%
%ex4:       
%           %load stimfile (e.g., stim01.mat) then run analyses
%           %defaults analyses are,'dataDis','signDistance','dataMeanAndStd'
%           analyses
%
%
%Mandatory inputs
%----------------
%
%                    subjects: {'sub01','sub02'} as many subjects as you want
%                         FAs: {'coh','Pstd','sample_dir'} for von Mises priors
%                              {'coh','priormodes','sample_dir'} for bimodal priors
%
%
%
%varargin
%--------
%
%    'dataPath','~/MydataDir': data directory
%'experiment','vonMisesPrior': experiment with von Mises priors
% 'experiment','bimodalPrior': experiment with bimodal priors
%            'dataMeanAndStd': draw data mean(std) in the different conditions
%                              with their sem.
%                   'dataRaw': draw the raw data sorted by predictors (conditions)
%                   'dataDis': draw data distribution per condition
%    'dataDis','signDistance': draw data distribution per condition
%                              against motion direction signed distance to prior
%     'dataDis','absDistance': draw data distribution per condition
%                              against motion direction absolute distance to prior
%        'overlay','subjects': overlay subjects data
%          'dataPastEffect01': test if data are biased toward the past trial
%                              (recency bias hypothesis)
%          'dataPastEffect02': test if data are biased toward the past trial
%                              (recency bias hypothesis, inspired by Whitney
%                              et al., 2014 with two main differences:
%                             - previous trial is never the prior mean in which case a bias
%                               toward the prior mean will be interpreted as a bias toward the
%                               previous trial
%
%                             - previous trial is never on the same side as the prior mean with
%                               respect to current for the same reason as above.
%                     'dataRT': draw reaction times per condition
%'dataRT','RTforPriorVsMotDir':
%'dataRT','RTPriorVsMotDirGrp': group estimated direction in two groups
%                               near prior vs near motion dir
%                   'fastrun': use databank saved in current directory instead
%                              or re-computing the whole databank again
%
%
%stuff to do
%-----------
%- add a function to get information about data (nb subj, sessions etc...)
%- add code lines to check if the FAs are correctly input(orthograph etc..)
%- allow to input only one or two FAs (not always three)
%- center mean and variability plot on the prior mean.
%
%notes: requires my code library
%
%
%
%examples
%unimodal (Exp00 - Gaussian prior std :[Inf,20], coh :[.12 0.35 1])
%------------------------------------------------------------------
%     analyses({'sub01','sub02','sub03','sub04','sub05','sub06','sub07'},...
%            {'coh','Pstd','sample_dir'},...
%            'dataPath','/Users/steeve/data/dataPsychophy/Exp00_Variance_PstdInf020_Pmean225_3coh_20dir_fdb_01_PwMate_AdjTrial',...
%            'experiment','vonMisesPrior',...
%            'dataMeanAndStd');
%
%unimodal final (Exp01 - Von Mises prior std :[10,20,40,80], coh :[.06 .12 .24])
%-------------------------------------------------------------------------------
%     analyses({'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub11','sub12'},...
%            {'coh','Pstd','sample_dir'},...
%            'dataPath','/Users/steeve/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%            'experiment','vonMisesPrior',...
%            'dataMeanAndStd');
%
%     analyses({'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},...
%            {'coh','Pstd','sample_dir'},...
%            'dataPath','/Users/steeve/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%            'experiment','vonMisesPrior',...
%            'dataPastEffect02');
%
%     analyses({'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},...
%            {'coh','Pstd','sample_dir'},...
%            'dataPath','/Users/steeve/Steeve_psychophy_data/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%            'experiment','vonMisesPrior',...
%            'dataPastEffect01');
%
%     analyses({'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},...
%            {'coh','Pstd','sample_dir'},...
%            'dataPath','/Users/steeve/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%            'experiment','vonMisesPrior',...
%            'dataDis','signDistance');
%
%     analyses({'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub11','sub12'},...
%            {'coh','Pstd','sample_dir'},...
%            'dataPath','/Users/steeve/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%            'experiment','vonMisesPrior',...
%            'dataDis','signDistance',...
%            'overlay','subjects');
%
%bimodal
%-------
% analyses({'sub03'},...
%     {'coh','priormodes','sample_dir'},...
%     'dataPath',...
%     '/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%     'experiment','bimodalPrior',...
%     'dataDis','signDistance')
%
% analyses({'sub03'},...
%     {'coh','priormodes','sample_dir'},...
%     'dataPath',...
%     '/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%     'experiment','bimodalPrior',...
%     'dataPastEffect02');
%
%     analyses({'sub01','sub02','sub03'},...
%            {'coh','priormodes','sample_dir'},...
%            'dataPath',...
%            '/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%            'experiment','bimodalPrior',...
%            'dataMeanAndStd');
%
%     analyses({'sub01','sub02','sub03'},...
%            {'coh','priormodes','sample_dir'},...
%            'dataPath',...
%            '/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%            'experiment','bimodalPrior',...
%            'dataDis','signDistance');
%
%     analyses({'sub01','sub02','sub03'},...
%            {'coh','priormodes','sample_dir'},...
%            'dataPath',...
%            '/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%            'experiment','bimodalPrior',...
%            'dataPastEffect01');
%
%     analyses({'sub01','sub02','sub03'},...
%            {'coh','priormodes','sample_dir'},...
%            'dataPath',...
%            '/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%            'experiment','bimodalPrior',...
%            'dataPastEffect02');
%



%Main analyses
%-------------
function [fig,GenStat] = analyses(subjects,FAs,varargin)

%------------------- Quick check data ------------------
%data are already loaded in the workspace
if nargin==0
    slPrintfStr('analyses',' Checking data...')   
    %get data       
    databank = SLMakedatabank;    
    %default factors   
    FAs = {'StimStrength','Pstd','FeatureSample'};
    %default analyses
    varargin = {'experiment','vonMisesPrior','dataDis','signDistance','dataMeanAndStd'};
    slPrintfStr('analyses',' Default factors are...{StimStrength,Pstd,FeatureSample};');   
    slPrintfStr('analyses',' Default analyses are...,"dataDis,signDistance,dataMeanAndStd"');   
else
    
%-------------------- subject-based analysis --------------
    
    %call for help
    if ieNotDefined('subjects')
        help analyses
        return
    end
    
    %get the analysis
    fprintf('Your analyses are: ')
    
    %set data path
    %datapath must always be input
    if sum(strcmp(varargin,'dataPath'))
        
        %case default run
        if ~sum(strcmp(varargin,'fastrun'))==1
            
            %Gather all available data in directory.
            dataPath = varargin{find(strcmp(varargin,'dataPath'))+1};
            cd(dataPath)
            fprintf('%s \n','(analyses) Now creating a databank...')
            databank = SLMakedatabank(subjects,varargin);
        else
            %case fast run.
            %Check if databank has already been created
            isCreatedDatabank = dir('*datbank.mat');
            
            %If not, create
            if isempty(isCreatedDatabank)==1
                fprintf('%s \n',['(analyses) Databank was not found in the'...
                    ' directory. Creating databank...'])
                dataPath = varargin{find(strcmp(varargin,'dataPath'))+1};
                cd(dataPath)
                databank = SLMakedatabank(subjects,varargin);
            else
                %If yes, load
                fprintf('%s \n',['(analyses) Databank was found in the'...
                    ' directory. Loading ...'])
                load datbank.mat
            end
        end
    else
        %case no datapath input
        fprintf('I couldn t find your path...please enter the path manually.')
        dataPath = uigetdir(cd,'Pick a project directory e.g., home/data/dataPsychophy/Exp0.../data');
        cd(dataPath)
        
        %create databank
        varargin{length(varargin)+1} = 'dataPath';
        varargin{length(varargin)+1} = dataPath;
        databank = SLMakedatabank(subjects,varargin);
    end
end

%outputs
GenStat = [];

%get the experimental conditions
fprintf('%s \n','(analyses) I am getting the factors...please wait')
[FA,est_dir] = initFAs(databank,FAs);


%check which analyses to run
%---------------------------
%look at data stats per condition
if sum(strcmp(varargin,'dataMeanAndStd'))
    [fig,GenStat] = plotMeanAndVarBehavior(databank,FA,est_dir,varargin);
end

%look at raw data per condition
if sum(strcmp(varargin,'dataRaw'))
    drawRaw(databank)
end

%look at data distribution per condition
if sum(strcmp(varargin,'dataDis'))
    fprintf('%s \n','(analyses) I am plotting the data distributions per condition...please wait')
    type = varargin(find(strcmp(varargin,'dataDis'))+1);
    dataDensity(databank,0.06,80,type,varargin)
end

%Look for "sequential effects 01"
if sum(strcmp(varargin,'dataPastEffect01'))
    dataPastEffect01(databank,FA,varargin)
end

%Look for "sequential effects 02"
if sum(strcmp(varargin,'dataPastEffect02'))
    dataPastEffect02(databank,FA,varargin)
end

%Look at reaction Times
if sum(strcmp(varargin,'dataRT'))
    
    %status
    fprintf('%s \n','(analyses) I am evaluating the reaction times...please wait')
    
    %analysis
    ReactionTime = cell2mat(databank.data(:,(strcmp(databank.nm,'ReactionTime'))==1));
    drawStat(ReactionTime,FA.g3.thisT,FA.g1.thisT,FA.g2.thisT);
    
    %label
    SLplotLabelAll('Reaction Times (s)','Motion directions (deg)')
    SLConventionUp
    fprintf('%s \n','All done.')
    
    %case check if reaction times differ when subjects estimate are near
    %prior versus when estimates are near motion directions. Mean reaction
    %times at each direction estimates.
    if sum(strcmp(varargin,'RTforPriorVsMotDir'))
        
        %ReactionTime
        %databank.estimatesDeg
        type = varargin(find(strcmp(varargin,'RTforPriorVsMotDir'))+1);
        drawReactionTimesPerCon(databank,0.06,80,type,varargin)
    end
    
    %Mean reaction times grouped based on distance to prior and motion direction.
    if sum(strcmp(varargin,'RTPriorVsMotDirGrp'))
        type = varargin(find(strcmp(varargin,'RTPriorVsMotDirGrp'))+1);
        drawRTPriorVsMotDirGrp(databank,type,varargin)
    end
end


fprintf('(analyses) Done ! \n')


%%factors (predictors)
function [FA,est_dir] = initFAs(databank,FAs)

%get data
%--------
%Get est data (e.g., estimated dirs)
est_data = cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));

%Set factors
%FA 1
FA.g1.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{1}))==1));
FA.g1.nm       =FAs{1};
%FA 2
FA.g2.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{2}))==1));
FA.g2.nm=FAs{2};
%FA 3
FA.g3.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{3}))==1));
FA.g3.nm       =FAs{3};


%CASE 1: Prior is unimodal
%-------------------------
%Get and order levels of each group
%group 1 (e.g., priors)
FA.g1.lvlsnm=unique(FA.g1.thisT); %names
FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend'); %order
FA.g1.lvlsnb=numel(FA.g1.lvlsnm);   %number
clear i
for i=1:FA.g1.lvlsnb
    index.g1.lvli(i)={find(FA.g1.thisT==FA.g1.lvlsnm(i))};
end
%group 2 (e.g., coherences)
FA.g2.lvlsnm=unique(FA.g2.thisT); %names
FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend'); %order
FA.g2.lvlsnb=numel(FA.g2.lvlsnm);   %number
for i=1:FA.g2.lvlsnb
    index.g2.lvli(i)={find(FA.g2.thisT==FA.g2.lvlsnm(i))};
end
%group 3 (e.g., displayed dirs)
FA.g3.lvlsnm=unique(FA.g3.thisT); %names
FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'ascend'); %order
FA.g3.lvlsnb=numel(FA.g3.lvlsnm);   %number
for i=1:FA.g3.lvlsnb
    index.g3.lvli(i)={find(FA.g3.thisT==FA.g3.lvlsnm(i))};
end


%CASE 2, Prior is bimodal
%------------------------
%factor 1 are the prior modes sorted based on how far they are from each
%other: the distance between the modes.

%if factor 1 is the prior
clear i
if strcmp(FA.g1.nm,'priormodes')
    
    %record the distance between the modes of each bimodal prior
    FA.g1.distance_mode=FA.g1.thisT(:,2)-FA.g1.thisT(:,1);
    FA.g1.lvlsnm=unique(FA.g1.distance_mode);
    FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend');
    FA.g1.lvlsnb=numel(FA.g1.lvlsnm);
    
    %find trials position of each prior condition
    for i=1:FA.g1.lvlsnb
        index.g1.lvli(i) = {find(FA.g1.distance_mode==FA.g1.lvlsnm(i))};
    end
end

%if factor 2 is the prior
clear i
if strcmp(FA.g2.nm,'priormodes')
    
    %record the distance between the modes of each bimodal prior
    FA.g2.distance_mode=FA.g2.thisT(:,2)-FA.g2.thisT(:,1);
    FA.g2.lvlsnm=unique(FA.g2.distance_mode);
    FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend');
    FA.g2.lvlsnb=numel(FA.g2.lvlsnm);
    
    %find trials position of each prior condition
    for i=1:FA.g2.lvlsnb
        index.g2.lvli(i)={find(FA.g2.distance_mode==FA.g2.lvlsnm(i))};
    end
    
end

%if factor 3 is the prior
clear i
if strcmp(FA.g3.nm,'priormodes')
    
    %record the distance between the modes of each bimodal prior
    FA.g3.distance_mode=FA.g3.thisT(:,2)-FA.g3.thisT(:,1);
    FA.g3.lvlsnm=unique(FA.g3.distance_mode);
    FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'descend');
    FA.g3.lvlsnb=numel(FA.g3.lvlsnm);
    
    %find trials position of each prior condition
    for i=1:FA.g3.lvlsnb
        index.g3.lvli(i)={find(FA.g3.distance_mode==FA.g3.lvlsnm(i))};
    end
end

%Calculate coordinates of average estimated directions for each condition
%organize data in following order: group 1(subplots) - group 2(colors) -
%group 3(x-axis). Each cell contains repetitions of a condition.
for k=1:FA.g1.lvlsnb
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            index.g1g2g3(j,i,k)={intersect( intersect( index.g1.lvli{k},index.g2.lvli{j} ),index.g3.lvli{i} )};
            
            %calculate statistics of estimated dirs (j:group 2, i:group 3, k:group 1)
            est_dir{j,i,k} = SLstatcircular(est_data(index.g1g2g3{j,i,k},:));
            
        end
    end
end

%%Test "Bayesian optimality" Hyp. with data average and variability
function [figs,GenStat] = plotMeanAndVarBehavior(databank,FA,est_dir,varargin)

%initialize figure
fig1.hdle=figure('color','w');
fig1.nm = ['myAnalysis','_mean'];
set(fig1.hdle, 'units', 'centimeters', 'pos', [0 0 30 10])

%graphic parameters
%one colormap
FA.g2.color={[0.5 0 0],...
    [1 0 0],...
    [1 0.5 0],...
    [0.65 0.65 0]};

%for fit
FA.g2.colorfit={[0.7 0.2 0],...
    [1 0.2 0],...
    [1 0.7 0],...
    [0.7 0.7 0]};

%(case bimodal prior) get unique prior modes
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
        priormodes=cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
        uniquePriormodes=unique(sort(priormodes,2), 'rows');
    end
end

%behavior average (mean)
%-----------------------
fprintf('%s \n','(plotMeanAndVarBehavior) Now plotting the mean of the estimated data...')
for k=1:FA.g1.lvlsnb
    subplot(1,FA.g1.lvlsnb,k)
    axis square
    hold all
    title([FA.g1.lvlsnm(k)])
    
    %(case Gaussian prior) indicate the priors' mean
    if sum(strcmp(varargin{1},'experiment'))
        if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
            %plot([FA.g3.lvlsnm(1) FA.g3.lvlsnm(end)],...
            %    [databank.data{1,7} databank.data{1,7}],...
            %    'b:',...
            %    'linewidth',2);
            plot([1 360],...
                [databank.data{1,7} databank.data{1,7}],...
                'b:',...
                'linewidth',1);
            hold on; plot([0 360],...
                [0 360],...
                'k:',...
                'linewidth',2);
        end
    end
    
    %(case bimodal prior) indicate the priors' modes
    if sum(strcmp(varargin{1},'experiment'))
        if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
            
            %plot the modes of each prior
            for z=1:size(uniquePriormodes,1)
                
                %get this prior's two modes and plot them
                thesemodes=repmat(uniquePriormodes(z,:),2,1);
                plot(repmat([FA.g3.lvlsnm(1),FA.g3.lvlsnm(end)]',1,2),thesemodes,...
                    ':',...
                    'color',FA.g2.color{z},...
                    'linewidth',1);
            end
        end
    end
    
    %indicate ideal estimations
    plot(FA.g3.lvlsnm,FA.g3.lvlsnm','k:','linewidth',2);
    
    %plot estimates stats for each experimental condition
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            
            %data
            fig1.datamean.deg(j,i,k)=est_dir{j,i,k}.deg.mean;
            fig1.datasem(j,i,k)=est_dir{j,i,k}.deg.sem;
            fig1.datamean.coord{j,i,k}=est_dir{j,i,k}.coord.mean;
            
            %factors g1,g2,g3
            fig1.dataInfog1{j,i,k}=FA.g1.lvlsnm(k);
            fig1.dataInfog2{j,i,k}=FA.g2.lvlsnm(j);
            fig1.dataInfog3{j,i,k}=FA.g3.lvlsnm(i);
            
            %motion directions in cartesian coordinates
            r=2.5;
            %fig1.displayed.coord{j,i,k}=polar2cartesian(fig1.dataInfog3{j,i,k},r);
            
            %Geometric method
            %Calculate distance between displayed and estimated dirs.
            %fig1.Disp_dataDistance{j,i,k}=SLvectors2signedAngle(fig1.datamean.coord{j,i,k}, fig1.displayed.coord{j,i,k} );
            
            %Select raw data
            fig1.datamean.deg(j,i,k)=fig1.datamean.deg(j,i,k);
        end
        
        
        %plot mean data
        scatter(FA.g3.lvlsnm,fig1.datamean.deg(j,:,k),...
            100,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',FA.g2.color{j},...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))))
        
        %plot sem
        SLerrorbar(FA.g3.lvlsnm,fig1.datamean.deg(j,:,k),...
            'yError',fig1.datasem(j,:,k),...
            'MarkerSize=1',['Color=',num2str(FA.g2.color{j})])
        
        %Linear fit
        linefit(FA.g3.lvlsnm',fig1.datamean.deg(j,:,k),1:10:360,FA.g2.colorfit{j});        
    end
    
    %graphics
    %--------
    set(gca,...
        'xtick', [65 225 385],'xticklabel',[-160 0 +160],...
        'ytick', [65 225 385],'yticklabel',[-160 0 +160],...
        'fontsize',13);
    xlim([0 360]);
    ylim([0 360]);
    
    %set labels
    if k==1
        ylabel('Estimated direction (deg)','fontsize',14);
    end
    xlabel('Displayed direction (deg)','fontsize',14);
    
    %remove dead space
    drawPublishAxis
end
%SLremoveDeadSpace



%behavior variability (std)
%-------------------------
%Draw
fprintf('%s \n','(plotMeanAndVarBehavior) Now plotting the std of the data...')
fig2.hdle=figure('color','w');
fig2.nm='MyAnalysis_std';
set(fig2.hdle,'units','centimeters','pos',[0 0 30 10])

%Set the number of bootstrap
numboot=50;
tic
maxPlot=[];
for k = 1 : FA.g1.lvlsnb
    subplot(1,FA.g1.lvlsnb,k)
    axis square
    hold all
    
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            
            %data
            fig2.datastd(j,i,k)=est_dir{j,i,k}.deg.std;
            fig2.datanb(j,i,k)=est_dir{j,i,k}.num;
            
            %bootstrap data (density of std)
            %initialize the data to bootstrap
            data2boot=[];
            fig2.databootstdSeries{j,i,k}=[];
            
            %collect the data to bootstrap (e.g. data at each displayed dir)
            fig2.datacoord(j,i,k)=est_dir{j,i,k}.coord;
            data2boot.data=fig2.datacoord(j,i,k).all;
            
            %create samples
            for jj=1:numboot
                
                %count the data
                data2boot.nb=size(data2boot.data,1);
                
                %if data are observed in x and number of sample >=2
                if isempty(data2boot.data)==0&&data2boot.nb>=2
                    
                    %bootstrap the data in x
                    %sample randomly a data value from x 'data2boot.nb' times
                    for iii=1:data2boot.nb
                        jit=ceil(rand*data2boot.nb);
                        data_boottmp{j,i,k}{jj}(iii,:)=data2boot.data(jit,:);
                        
                        %store for fig2
                        fig2.databoot{j,i,k}{jj}(iii,:)=data_boottmp{j,i,k}{jj}(iii,:) ;
                        fig2.databootnb{j,i,k}{jj}(iii,:)=data2boot.nb;
                    end
                    
                    %calculate the std for each bootstrap
                    fig2.databootstd{j,i,k}{jj}=SLstatcircular(fig2.databoot{j,i,k}{jj});
                    
                    %collect the std in a vector
                    fig2.databootstdSeries{j,i,k}=[fig2.databootstdSeries{j,i,k} fig2.databootstd{j,i,k}{jj}.deg.std];
                    
                    %if there are no data in x
                else
                    %fill in a nan
                    data_boottmp{j,i,k}{jj}=nan;
                    
                    %store for fig2
                    fig2.databoot{j,i,k}{jj}=data_boottmp{j,i,k}{jj} ;
                    
                    %there are no data
                    fig2.databootnb{j,i,k}{jj}=0 ;
                    
                    %nan the std
                    fig2.databootstd{j,i,k}{jj}=nan;
                    
                    %nan the vector of std
                    fig2.databootstdSeries{j,i,k}=nan;
                end
            end
            
            %calculate the std over bootstrapped samples
            fig2.databootMeanOfstd{j,i,k}=nanmean(fig2.databootstdSeries{j,i,k});
            fig2.databootStdOfstd{j,i,k}=std(fig2.databootstdSeries{j,i,k});
            
            %factors g1,g2,g3
            fig2.dataInfog1{j,i,k}=FA.g1.lvlsnm(k);
            fig2.dataInfog2{j,i,k}=FA.g2.lvlsnm(j);
            fig2.dataInfog3{j,i,k}=FA.g3.lvlsnm(i);
        end
        
        %draw variability
        myerrorbar(FA.g3.lvlsnm,[fig2.databootMeanOfstd{j,:,k}], 'yError',...
            [fig2.databootStdOfstd{j,:,k}],...
            'Symbol=.',...
            'LineWidth=1',...
            ['Color=',num2str(FA.g2.color{j})],...
            ['MarkerFacecolor=', num2str(FA.g2.color{j})],...
            'MarkerSize=1',...
            'MarkerEdgeColor=w');
        
        scatter(FA.g3.lvlsnm,[fig2.databootMeanOfstd{j,:,k}],100,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',FA.g2.color{j},...
            'displayname',strcat(FA.g2.nm,':',num2str(FA.g2.lvlsnm(j))))
        
        h(k)=gca;
        maxPlot=[maxPlot max(sum([[fig2.databootMeanOfstd{j,:,k}]; [fig2.databootStdOfstd{j,:,k}]]))];
    end
    
    %graphics
    %--------
    %x-axis
    set(gca,'xtick',0:45:360,'xticklabel',0:45:360,...
        'fontsize',13)
    set(gca,'ytick',0:10:100,'yticklabel',0:10:100,'fontsize',13);
    xlim([0 360]);
    ylim([0 100]);
    
    %set labels of axes
    if k==1
        ylabel('Std of estimated direction (degrees)','fontsize',11);
    end
    xlabel({'Displayed direction (degrees)','(degrees)'},'fontsize',11);
    
    %make clean
    drawPublishAxis
end
%SLremoveDeadSpace

% for k=1:FA.g1.lvlsnb
%     ylim(h(k),[0 max(maxPlot)])
%     set(h(k),'ytick',1:10:max(maxPlot),'yticklabel',1:10:max(maxPlot))
% end


%Draw raw data
% fprintf('%s \n','(drawRawEstDir) Now plotting individual estimates in polar...')
% fig4=drawRawEstDir(fig2);

%backup figures informations
figs={fig1, fig2};


%%Stat of ANCOVA: effect of FA1(e.g., coherence) on the slope of data=f(FA3)).
%%----------------------------------------------------------------------
%%Set x: (e.g.,FA3 - displayed dirs)
%x=[];
%x=reshape(repmat(FA.g3.lvlsnm , FA.g2.lvlsnb*FA.g1.lvlsnb,1), [],1);
%
%%Set y: (data; e.g., estimated dirs)
%%Loop over the levels of Factor 1
%y=[];
%for k=1:FA.g1.lvlsnb
%    ytmp=[];
%    ytmp=reshape(fig1.datamean.deg(:,:,k)',[],1);
%    y=[y; ytmp];
%end
%
%%Set g the groups (e.g., FA1 - coherences)
%g=[];
%g=reshape(repmat(1: FA.g1.lvlsnb, FA.g3.lvlsnb*FA.g2.lvlsnb,1), [],1);
%
%%Run the Ancova
%Ancova_effect_of_FA1=statAncova(x, y, g);
%display('Your Ancova is being performed...')
%
%%%Report the statistics of the Ancova
%%%Degrees of freedom
%%%Within groups  (e.g.,=f(nb of subjects))
%%dfwgT1  =GenStat{2,1}{2,2}{5,2};
%%%Between groups  (e.g.,=f(nb of group))
%%dfbgT1  =GenStat{2,1}{2,2}{4,2};
%%%F-value for the interaction effect
%%FvalueT1=GenStat{2,1}{2,2}{4,5};
%%%p-value of the tested interaction effect
%%pvalueT1=GenStat{2,1}{2,2}{4,6};
%%GenStat2ReportT1={dfbgT1, dfwgT1, FvalueT1, pvalueT1};
%%display(GenStat2ReportT1)

%
%----------------------------------------------------------------------
%Stat of ANCOVA: effect of FA1(e.g., prior's width) on the slope of data=f(FA3)).
%----------------------------------------------------------------------
%Set x: (e.g.,FA3 - displayed dirs)
x=[];
x=reshape(repmat(FA.g3.lvlsnm , FA.g2.lvlsnb*FA.g1.lvlsnb,1), [],1);

%Set y: (data; e.g., estimated dirs)
%Loop over the levels of Factor 1
y=[];
for k=1:FA.g1.lvlsnb
    ytmp=[];
    ytmp=reshape(fig1.datamean.deg(:,:,k)',[],1);
    y=[y; ytmp];
end

%Set g the groups (e.g., FA1 - coherences)
g=[];
g=reshape(repmat(1: FA.g2.lvlsnb, FA.g3.lvlsnb*FA.g1.lvlsnb,1), [],1);

%check enough groups
if length(unique(g))<2
    GenStat = [];
    sprintf('(plotMeanAndVarBehavior) Not enough groups for ANCOVA')
    return
end

%Run the Ancova
Ancova_effect_of_FA2 = statAncova(x, y, g);
fprintf('%s \n','(analyses) running Ancova...')

%%Report the statistics of the Ancova
%%Degrees of freedom
%%Within groups  (e.g.,=f(nb of subjects))
%dfwgT1  =GenStat{2,1}{2,2}{5,2};
%%Between groups  (e.g.,=f(nb of group))
%dfbgT1  =GenStat{2,1}{2,2}{4,2};
%%F-value for the interaction effect
%FvalueT1=GenStat{2,1}{2,2}{4,5};
%%p-value of the tested interaction effect
%pvalueT1=GenStat{2,1}{2,2}{4,6};
%GenStat2ReportT1={dfbgT1, dfwgT1, FvalueT1, pvalueT1};
%display(GenStat2ReportT1)


%%----------------------------------------------------------------------
%Stat T1: ANCOVA (effect of FA2 on the slope of data=f(FA3) compared to 1.
%----------------------------------------------------------------------
%e.g., FA 1 is coherence
%e.g., FA 2 is prior's strength
%e.g., FA 3 is displayed dir

%Loop over each condition FA1-FA2
%Levels of FA 1
for k=1:FA.g1.lvlsnb
    %Levels of FA 2
    for j=1:FA.g2.lvlsnb
        %Set x: displayed dirs (levels of FA3)
        x=[];
        x=reshape(repmat(FA.g3.lvlsnm, 2, 1), [],1);
        %Set y: data (e,g., estimated dirs and dirs with slope 1).
        y=[];
        y=reshape([fig1.datamean.deg(j,:,k)', FA.g3.lvlsnm], [],1);
        %Set the groups (g)
        g=[];
        g=reshape(repmat([1,2], FA.g3.lvlsnb,1), [],1);
        %    %Title rows
        %    GenStatAgainst1{j+1,  1}=[FA.g2.nm num2str(FA.g2.lvlsnm(j))];
        %    %title columns
        %    GenStatAgainst1{  1,k+1}=[FA.g1.nm num2str(FA.g1.lvlsnm(k))];
        %run analysis
        Ancova_effect_of_FA2vs1{j+1,k+1}=statAncova(x,y,g);
    end
end

%Stat to report (write sthg here)


%%---------------------------------------------------------------------
%Stat: ANCOVA - compare the slopes of FA 2's levels for each FA 1.
%---------------------------------------------------------------------
%e.g., FA 1 is coherence
%e.g., FA 2 is prior's strength
%e.g., FA 3 is displayed dir
%Loop over the levels of FA1 (e.g., coherences)
for k=1:FA.g1.lvlsnb
    %Set x: FA3 - e.g., displayed dirs (levels of FA 3)
    x = [];
    x = reshape(repmat(FA.g3.lvlsnm, FA.g2.lvlsnb,1),[],1);
    
    %Set y: data e.g., estimated dirs
    y = [];
    y = reshape(fig1.datamean.deg(:,:,k)', [], 1);
    
    %Set g: groups (levels of FA2, e.g., priors)
    g = [];
    g = reshape(repmat(FA.g2.lvlsnm', FA.g3.lvlsnb,1), [],1);
    
    %title columns
    GenStatBetwLv{1,k} = [FA.g1.nm num2str(FA.g1.lvlsnm(k))];
    
    %run analysis
    GenStatBetwLv{2,k} = statAncova(x,y,g);
end

%Stat to report(write sthg here)


%----------------------------------------------------------------------
%Store the statistics
%----------------------------------------------------------------------
%Title columns
GenStat=[ {'Ancova_effect_of_FA2'},...
    {'Ancova_effect_of_FA2vs1'},...
    {'Ancova_lvl1_vs_lvl2_FA1'}];
%Data
GenStat(2,:)=[ {Ancova_effect_of_FA2},...
    {Ancova_effect_of_FA2vs1},...
    {GenStatBetwLv}];
fprintf ('%s \n','(plotMeanAndVarBehavior) The results of the Ancova have been stored')

%Report the statistics of the Ancova
%Degrees of freedom
%Within groups  (e.g.,=f(nb of subjects))
dfwgT1  =GenStat{2,1}{2,2}{5,2};
%Between g roups  (e.g.,=f(nb of group))
dfbgT1  =GenStat{2,1}{2,2}{4,2};
%F-value for the interaction effect
FvalueT1=GenStat{2,1}{2,2}{4,5};
%p-value of the tested interaction effect
pvalueT1=GenStat{2,1}{2,2}{4,6};
GenStatT1={dfbgT1, dfwgT1, FvalueT1, pvalueT1};
display(GenStatT1)


%Report the statistics of the Ancova
%Degrees of freedom
%Within groups  (e.g.,=f(nb of subjects))
dfwgT2  =GenStat{2,2}{2,2}{2,2}{5,2};
%Between groups  (e.g.,=f(nb of group))
dfbgT2  =GenStat{2,2}{2,2}{2,2}{4,2};
%F-value for the interaction effect
FvalueT2=GenStat{2,2}{2,2}{2,2}{4,5};
%p-value of the tested interaction effect
pvalueT2=GenStat{2,2}{2,2}{2,2}{4,6};
GenStatT2={dfbgT2, dfwgT2, FvalueT2, pvalueT2};
display(GenStatT2)


%Report the statistics of the Ancova
%Degrees of freedom
%Within groups  (e.g.,=f(nb of subjects))
dfwgT3  =GenStat{2,3}{2,1}{2,2}{5,2};
%Between groups  (e.g.,=f(nb of group))
dfbgT3  =GenStat{2,3}{2,1}{2,2}{4,2};
%F-value for the interaction effect
FvalueT3=GenStat{2,3}{2,1}{2,2}{4,5};
%p-value of the tested interaction effect
pvalueT3=GenStat{2,3}{2,1}{2,2}{4,6};
GenStatT3={dfbgT3, dfwgT3, FvalueT3, pvalueT3};
display(GenStatT3)



%Backup figure
autobackup(fig1.hdle,fig1.nm,'.fig');

%[fig]=graphStat(FA, est_dir,fig);
%Graph the quantitative effect of the prior's strength
if FA.g1.lvlsnb ==2
    [fig3]=graphPriorEffect(fig1,FA);
end

%%%Test "Bayesian optimality" Hyp. with data density
function dataDensity(databank,StimStrengthi,pstdi,type,varargin)

%%get data
dir.es = cell2mat(databank.data(:,strcmp(databank.nm,'estimatedFeature')==1));

%experimental factors
%prior
prior.m = unique(cell2mat(databank.data(:,strcmp(databank.nm,'priormean')==1)));
pstd = cell2mat(databank.data(:,strcmp(databank.nm,'Pstd')==1));

%motion directions
dir.di = cell2mat(databank.data(:,strcmp(databank.nm,'FeatureSample')==1));
x = [0 5:10:360];

%coherence
StimStrength = cell2mat(databank.data(:,strcmp(databank.nm,'StimStrength')==1));

%(case gaussian prior)
%---------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
        
        %draw data sorted by displayed directions with same abs(distance) to prior mean
        if strcmp(type,'absDistance')
            
            %graphics
            %densities
            %cmp=colormap('jet');
            
            %color leftward density
            colors{1}=[.5 .5 .5];
            
            %color rightward density
            colors{2}=[.9 0.8 0.5];
            
            %markers
            colormk(1,:)=[0 0 0];
            colormk(2,:)=[.7 0.6 0.5];
            
            %Organize data for a clear plot
            %The data for diplayed direction that are equidistant to the prior are
            %displayed on the same subplot.
            %sort data
            C.StimStrengthPstd=StimStrength==StimStrengthi&pstd==pstdi;
            
            %range of displayed directions
            dir.diStimStrengthPrior=dir.di(C.StimStrengthPstd);
            dir.di=unique(dir.diStimStrengthPrior);
            dir.dicount=numel(dir.di);
            
            %get distances to prior
            dist2Pr.sign=dir.di-225;
            dist2Pr.abs=abs(dist2Pr.sign);
            dist2Pr.uniq=unique(dist2Pr.abs);
            dist2Pr.uniqCount=numel(dist2Pr.uniq);
            
            %Initialize figure
            scrsz=get(0,'ScreenSize');
            figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)])
            hold all
            subp=nan(dist2Pr.uniqCount,1);
            for i=1:dist2Pr.uniqCount
                subp(i)=subplot(dist2Pr.uniqCount,1,i);
            end
            dir.dist2Pr=dir.diStimStrengthPrior-225;
            
            %estimated directions
            dir.esStimStrengthPrior=dir.es(C.StimStrengthPstd);
            
            %Look at each displayed direction
            dirleftandright=nan(2,1);
            for i=1:dist2Pr.uniqCount
                
                %draw
                subplot(subp(i))
                
                %directions equidistant to prior are displayed together
                %we find the direction to the left of the prior, then to the one to the
                %right of the prior and we plot them.
                eqpos=find(dist2Pr.abs==dist2Pr.uniq(i));
                pmax=nan(numel(eqpos));
                for j=1:numel(eqpos)
                    
                    %sort data
                    dir.dispos=dir.dist2Pr==dist2Pr.sign(eqpos(j));
                    C.dispStimStrengthPrior=dir.dispos;
                    dirleftandright(j)=unique(dir.diStimStrengthPrior(C.dispStimStrengthPrior));
                    
                    %if this condition exists
                    if sum(C.dispStimStrengthPrior)~=0;
                        
                        %calculate pdf
                        [p,xpdf]=makePdf(x,dir.esStimStrengthPrior(C.dispStimStrengthPrior),'raw');
                        pmax(j)=max(p);
                        
                        %show sample size
                        dir.esCount=numel(dir.esStimStrengthPrior(C.dispStimStrengthPrior));
                        fprintf('%s %12g \n', '(analyses) number of trial is', dir.esCount)
                        
                        %draw
                        hold all
                        area(xpdf,p,'facecolor',colors{j},'edgecolor','none');
                        
                        %motion directions
                        %to the left then to the right
                        if ~isnan(dirleftandright(j))
                            plot([dirleftandright(j) dirleftandright(j)],[0 0.5],'-','color',colormk(j,:),'linewidth',3);
                        end
                        
                        %prior mean
                        plot([prior.m prior.m],[0 0.5],'-','color',[0.3 0.6 0.9],'linewidth',3);
                        axis off
                    end
                end
                
                %put on top
                %uistack(h1,'top')
                %uistack(h2,'top')
                
                %label
                set(gca,'xtick',[],'xticklabel',[])
                
                %title the
                if i==1
                    title([num2str(StimStrengthi*100),'% Stimulus strength and prior std = ',num2str(pstdi)],'fontsize',11)
                end
                if i==dist2Pr.uniqCount
                    ylabel(subp(i),'Probability','fontsize',11)
                    xlabel({'Estimated directions','(degrees)'},'fontsize',11)
                    set(gca,'xtick',[65 225 385],'xticklabel',[-160 0 +160])
                end
                xlim([0 400])
                ylim([0 max(max(pmax))])
                %drawPublishAxis
            end
        end
        
        %draw data sorted by signed distance of displayed direction to prior mean
        if strcmp(type,'signDistance')
            
            %plot color
            %colors=[.9 0.8 0.5];
            colors=[1 0 0];
            
            %set of experimental priors
            Thepriors=unique(pstd);
            
            %set of stimulus strengths
            TheStimStrengths=unique(StimStrength);
            
            %motion direction as signed linear distance to prior mean
            dir.dilinDisttoPrior=dir.di-prior.m;
            
            %set of motion directions as signed linear distance to prior mean
            ThedirlinDisttoPrior=unique(dir.dilinDisttoPrior);
            
            %look at how data distribution changes as prior strength increases
            %get number of axes
            axesPos=1:1:numel(ThedirlinDisttoPrior)*numel(Thepriors);
            dirpos=[];
            
            %calculate position of plot for each motion direction
            for i=1:numel(Thepriors)
                dirpos=[dirpos; i:numel(Thepriors):numel(axesPos)];
            end
            
            %look at this stimulus strength
            for j = 1 : numel(TheStimStrengths)
                thisStimStrength = TheStimStrengths(j);
                
                %open figure
                scrsz = get(0,'ScreenSize');
                width = 7;
                figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
                
                %look at this prior
                for i = 1 : numel(Thepriors)
                    thisprior = Thepriors(i);
                    
                    %Look at this motion direction
                    for k = 1 : numel(ThedirlinDisttoPrior)
                        
                        %get this motion direction
                        %in term of distance to the prior mean
                        thisDirlinDisttoPrior = ThedirlinDisttoPrior(k);
                        
                        %subplots
                        subplot(numel(ThedirlinDisttoPrior),numel(Thepriors),dirpos(i,k));
                        
                        %get data for this condition
                        dirtoHist = dir.es(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        
                        %if this condition exists
                        if ~isempty(dirtoHist)
                            
                            %calculate pdf
                            [p,xpdf] = makePdf(x,dirtoHist,'raw');
                            
                            %show sample size
                            fprintf('%s %12g \n', '(analyses) number of trial is',numel(dirtoHist))
                            
                            %draw
                            hold all
                            area(xpdf,p,'facecolor',colors,'edgecolor','none');
                            
                            %graphics
                            xlim([0 360])
                            ylim([0 max(p)])
                            axis tight
                            set(gca,'ytick',[0 fix(max(p)*100)/100],'ytickLabel',[0 fix(max(p)*100)/100])
                            
                            %indicate motion direction on the plot
                            %get motion direction
                            thisDir=thisDirlinDisttoPrior+prior.m;
                            plot([thisDir thisDir],[0 0.5*max(p)],'-','color',[0 0 0],...
                                'linewidth',3)
                            
                            %prior mean
                            plot([prior.m prior.m],[0 0.5*max(p)],'-','color',...
                                [0.3 0.6 0.9],'linewidth',3);
                            
                            %remove axis for clarity
                            axis off
                            set(gca,'ytick',0,...
                                'xticklabel',0)
                        else
                            axis off
                        end
                        box off
                        
                        %graphics
                        %--------
                        %title top axes
                        if k==1
                            title(['(',num2str(thisStimStrength*100),'% Stimulus strength, std of prior=',num2str(thisprior),')'],'fontsize',11)
                        end
                        %no x-tick
                        set(gca,'xtick',[],'xtickLabel',[],'fontsize',10)
                        
                        %x-axis only for the bottom axis
                        if k==numel(ThedirlinDisttoPrior) && i==1
                            ylabel('Probability','fontsize',11)
                        end
                        
                        %ylabel for one of the axis only
                        if k==numel(ThedirlinDisttoPrior)
                            axis on
                            xlabel('Estimated directions (degrees)','fontsize',11)
                            xlim([0 360])
                            set(gca,'xtick',[65 225 385],'xtickLabel',[65 225 25])
                        end
                    end
                end
            end
            
            %case we want to overlay subjects data
            if strcmp(varargin{1}(find(strcmp(varargin{1},'overlay'))+1),'subjects')
                
                %subjects
                sub=unique(databank.subjects);
                numSub=length(sub);
                
                %plot color
                colors=linspecer(numSub);
                
                %priors and stimulus strengths
                Thepriors=unique(databank.Pstd);
                TheStimStrengths=unique(databank.StimStrengtherence);
                
                %motion direction (signed linear distance to prior mean)
                dir.dilinDisttoPrior=dir.di-prior.m;
                
                %motion directions (signed linear distance to prior mean)
                ThedirlinDisttoPrior=unique(dir.dilinDisttoPrior);
                
                %look at how data distribution changes as prior strength
                %increases.
                %calculate position of plot for each motion direction
                axesPos=1:1:numel(ThedirlinDisttoPrior)*numel(Thepriors);
                dirpos=[];
                for i=1:numel(Thepriors)
                    dirpos=[dirpos; i:numel(Thepriors):numel(axesPos)];
                end
                
                %look at this stimulus strength
                for j=1:numel(TheStimStrengths)
                    
                    %figure
                    scrsz=get(0,'ScreenSize');
                    width=7;
                    figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
                    
                    %this stimulus strength
                    thisStimStrength = TheStimStrengths(j);
                    
                    %look at this prior
                    for i=1:numel(Thepriors)
                        thisprior=Thepriors(i);
                        
                        %Look at this motion direction
                        for k=1:numel(ThedirlinDisttoPrior)
                            
                            %get this motion direction
                            %in term of distance to the prior mean
                            thisDirlinDisttoPrior=ThedirlinDisttoPrior(k);
                            
                            %look this subject
                            for isub=1:length(sub)
                                
                                %subplots
                                subplot(numel(ThedirlinDisttoPrior),...
                                    numel(Thepriors),dirpos(i,k));
                                
                                %colors
                                set(gca,'NextPlot','replacechildren',...
                                    'ColorOrder',colors);
                                
                                %get data for this condition
                                dirtoHist=dir.es(pstd==thisprior & ...
                                    StimStrength==thisStimStrength & ...
                                    dir.dilinDisttoPrior==thisDirlinDisttoPrior ...
                                    & databank.subjects == isub);
                                
                                %if this condition exists
                                if ~isempty(dirtoHist)
                                    
                                    %axis
                                    axis on
                                    
                                    %calculate pdf
                                    [p,xpdf]=makePdf(x,dirtoHist,'raw');
                                    
                                    %show sample size
                                    fprintf('%s %12g \n',...
                                        '(analyses) number of trial is',...
                                        numel(dirtoHist))
                                    
                                    %draw
                                    hold all
                                    %area(xpdf,p,'facecolor',colors(isub,:),'edgecolor','none');
                                    plot(xpdf,p,'color',colors(isub,:),...
                                        'linesmoothing','on',...
                                        'linewidth',2);
                                    
                                    %graphics
                                    xlim([0 360])
                                    ylim([0 max(p)])
                                    axis tight
                                    set(gca,'ytick',[0 fix(max(p)*100)/100],...
                                        'ytickLabel',[0 fix(max(p)*100)/100])
                                    
                                    %indicate motion direction on the plot
                                    %get motion direction
                                    thisDir=thisDirlinDisttoPrior+prior.m;
                                    plot([thisDir thisDir],[0 0.5*max(p)],...
                                        '-','color',[0 0 0],...
                                        'linewidth',3)
                                    
                                    %prior mean
                                    plot([prior.m prior.m],[0 0.5*max(p)],...
                                        '-','color',...
                                        [0.3 0.6 0.9],'linewidth',3);
                                    
                                    %graphics
                                    %--------
                                    %no x-tick
                                    set(gca,'xtick',[],'xtickLabel',[],...
                                        'fontsize',10)
                                    
                                    %x-axis only for the bottom axis
                                    if k==round(numel(ThedirlinDisttoPrior)/2) && i==1
                                        ylabel('Probability','fontsize',11)
                                    end
                                    
                                    %ylabel for one of the axis only
                                    if k==numel(ThedirlinDisttoPrior)
                                        xlim([0 360])
                                        set(gca,'xtick',[65 225 385],...
                                            'xticklabel',[65 225 25])
                                        xlabel('Estimated directions (degrees)',...
                                            'fontsize',11)
                                    end
                                    %box off
                                    set(gca,'ytick',0,...
                                        'yticklabel',0)
                                else
                                    if k~=numel(ThedirlinDisttoPrior)
                                        %axis off
                                    else
                                        xlim([0 360])
                                        set(gca,'xtick',[65 225 385],...
                                            'xticklabel',[65 225 25])
                                        set(gca,'ytick',0,...
                                            'yticklabel',0)
                                        xlabel('Estimated directions (degrees)',...
                                            'fontsize',11)
                                    end
                                end
                                
                                %title top axes
                                if k==1
                                    title(['(',num2str(thisStimStrength*100),...
                                        '% stimulus strength, std of prior=',...
                                        num2str(thisprior),')'],...
                                        'fontsize',11)
                                end
                            end
                        end
                    end
                    
                    %remove dead space
                    %SLremoveDeadSpace;
                end
            end
        end
    end
end



%(case bimodal prior)
%--------------------
%The two modes are indicated for each displayed direction.
%get unique prior modes
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
        
        %get trial prior modes
        priormodes=cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
        
        %get experimental priors
        Thepriors=unique(sort(priormodes,2), 'rows');
        
        %draw data sorted by signed distance of displayed direction to prior mean
        %for visualization.
        if strcmp(type,'signDistance')
            
            %plot color
            %colors=[.9 0.8 0.5];
            colors=[1 0 0];
            
            
            %set of stimulus strength
            TheStimStrengths=unique(StimStrength);
            
            %set of motion directions
            Thedir=unique(dir.di);
            
            %motion direction as signed linear distance to prior mean
            dir.dilinDisttoPrior=dir.di-prior.m;
            
            %set of motion directions as signed linear distance to prior mean
            ThedirlinDisttoPrior=unique(dir.dilinDisttoPrior);
            
            %look at how data distribution changes as prior strength increases
            %get number of axes
            axesPos=1:1:numel(ThedirlinDisttoPrior)*numel(Thepriors);
            dirpos=[];
            
            %calculate position of the plot for each motion direction
            numPriors=size(Thepriors,1);
            dirpos=[];
            for i=1:numPriors
                dirpos=[dirpos; i:numPriors:numel(axesPos)];
            end
            
            %look at this stimulus strength
            for j=1:numel(TheStimStrengths)
                thisStimStrength=TheStimStrengths(j);
                
                %open figure
                scrsz=get(0,'ScreenSize');
                width=7;
                figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
                
                %look at this prior
                for i=1:numPriors
                    thisprior=Thepriors(i,:);
                    
                    %Look at this motion direction
                    for k=1:numel(ThedirlinDisttoPrior)
                        
                        %get this motion direction
                        %in term of distance to the prior mean
                        thisDirlinDisttoPrior=ThedirlinDisttoPrior(k);
                        
                        %subplots
                        subplot(numel(ThedirlinDisttoPrior),numPriors,dirpos(i,k));
                        
                        %get data for this condition
                        dirtoHist=dir.es(priormodes(:,1)==thisprior(1) & priormodes(:,2)==thisprior(2) & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        
                        %if data for this condition exists
                        if ~isempty(dirtoHist)
                            
                            %calculate pdf
                            [p,xpdf]=makePdf(x,dirtoHist,'raw');
                            
                            %show number of sample in this condition
                            fprintf('%s \n', ['(dataDensity) Prior:[',...
                                num2str(thisprior),'] - ',...
                                'Stim strength:[',num2str(thisStimStrength),'] - ',...
                                num2str(numel(dirtoHist)),' trials'])
                            
                            %draw
                            hold all
                            area(xpdf,p,'facecolor',colors,'edgecolor','none');
                            
                            %graphics
                            xlim([0 360])
                            ylim([0 max(p)])
                            axis tight
                            set(gca,'ytick',[0 fix(max(p)*100)/100],...
                                'ytickLabel',[0 fix(max(p)*100)/100])
                            
                            %indicate motion direction on the plot
                            %get motion direction
                            thisDir=thisDirlinDisttoPrior+prior.m;
                            plot([thisDir thisDir],[0 0.5*max(p)],...
                                '-','color',[0 0 0],...
                                'linewidth',3)
                            
                            %prior modes
                            plot([thisprior(1) thisprior(1)],...
                                [0 0.5*max(p)],'-','color',...
                                [0.3 0.6 0.9],'linewidth',3);
                            plot([thisprior(2) thisprior(2)],...
                                [0 0.5*max(p)],'-','color',...
                                [0.3 0.6 0.9],'linewidth',3);
                        else
                            plot(x,nan)
                            xlim([0 360])
                        end
                        box off
                        
                        %graphics
                        %--------
                        %title top axes
                        if k==1
                            title(['(',num2str(thisStimStrength*100),...
                                '% stim strength, std of prior=',...
                                num2str(thisprior),')'],'fontsize',11)
                        end
                        %no x-tick
                        set(gca,'xtick',[65 225 385],'xtickLabel',[],'fontsize',10)
                        
                        %x-axis only for the bottom axis
                        if k==numel(ThedirlinDisttoPrior) && i==1
                            ylabel('Probability','fontsize',11)
                        end
                        
                        %ylabel for one of the axis only
                        if k==numel(ThedirlinDisttoPrior)
                            xlabel('Estimated directions (degrees)','fontsize',11)
                            set(gca,'xtick',[65 225 385],'xticklabel',[65 225 25])
                        end
                    end
                end
            end
        end
    end
end

%%%Reaction times per estimates
function drawReactionTimesPerCon(databank,StimStrengthi,pstdi,type,varargin)

%%get data

%experimental factors
%prior
prior.m = unique(cell2mat(databank.data(:,strcmp(databank.nm,'priormean')==1)));
pstd=cell2mat(databank.data(:,strcmp(databank.nm,'Pstd')==1));

%motion directions
dir.di=cell2mat(databank.data(:,strcmp(databank.nm,'sample_dir')==1));
x=5:10:360;

%stimulus strength
StimStrength=cell2mat(databank.data(:,strcmp(databank.nm,'StimStrength')==1));

%(case von Mises prior)
%---------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
        
        %draw data sorted by displayed directions with same abs(distance) to prior mean
        if strcmp(type,'absDistance')
            
            %graphics
            %densities
            %cmp=colormap('jet');
            
            %color leftward density
            colors{1}=[.5 .5 .5];
            
            %color rightward density
            colors{2}=[.9 0.8 0.5];
            
            %markers
            colormk(1,:)=[0 0 0];
            colormk(2,:)=[.7 0.6 0.5];
            
            %Organize data for a clear plot
            %The data for diplayed direction that are equidistant to the prior are
            %displayed on the same subplot.
            %sort data
            C.StimStrengthPstd = StimStrength==StimStrengthi&pstd==pstdi;
            
            %range of displayed directions
            dir.diStimStrengthPrior = dir.di(C.StimStrengthPstd);
            dir.di = unique(dir.diStimStrengthPrior);
            dir.dicount = numel(dir.di);
            
            %get distances to prior
            dist2Pr.sign=dir.di-225;
            dist2Pr.abs=abs(dist2Pr.sign);
            dist2Pr.uniq=unique(dist2Pr.abs);
            dist2Pr.uniqCount=numel(dist2Pr.uniq);
            
            %Initialize figure
            scrsz = get(0,'ScreenSize');
            figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)])
            hold all
            subp = nan(dist2Pr.uniqCount,1);
            for i = 1 : dist2Pr.uniqCount
                subp(i) = subplot(dist2Pr.uniqCount,1,i);
            end
            dir.dist2Pr = dir.diStimStrengthPrior - 225;
            
            %estimated directions and reaction times in this condition
            dir.esStimStrengthPrior = databank.estimatesDeg(C.StimStrengthPstd);
            dir.RTStimStrengthPrior = databank.reactionTime(C.StimStrengthPstd);
            
            %Look at each displayed direction
            dirleftandright = nan(2,1);
            for i = 1 : dist2Pr.uniqCount
                
                %draw
                subplot(subp(i))
                
                %directions equidistant to prior are displayed together
                %we find the direction to the left of the prior, then to the one to the
                %right of the prior and we plot them.
                eqpos = find(dist2Pr.abs==dist2Pr.uniq(i));
                pmax = nan(numel(eqpos));
                for j = 1 : numel(eqpos)
                    
                    %sort data per displayed direction distance to prior
                    %dir.dispos = dir.dist2Pr==dist2Pr.sign(eqpos(j));
                    C.dispStimStrengthPrior = dir.dist2Pr==dist2Pr.sign(eqpos(j));
                    
                    %if this condition exists
                    if sum(C.dispStimStrengthPrior)~=0
                        
                        %calculate pdf
                        %[p,xpdf]=makePdf(x,dir.esStimStrengthPrior(C.dispStimStrengthPrior),'raw');
                        stats = SLmakeStat(dir.RTStimStrengthPrior(C.dispStimStrengthPrior),dir.esStimStrengthPrior(C.dispStimStrengthPrior));
                        pmax(j) = max(stats.mean);
                        
                        %show sample size
                        RTCount = sum(~isnan(dir.RTStimStrengthPrior(C.dispStimStrengthPrior)));
                        fprintf('%s %12g \n', '(analyses) number of trial is',RTCount)
                        
                        %----
                        %draw
                        %-----
                        hold all
                        %area(stats.conditions(:,1),stats.mean,'facecolor',colors{j},'edgecolor','none');
                        y = nan(360,1);
                        ystd = nan(360,1);
                        y(stats.conditions(:,1)) = stats.mean;
                        ystd(stats.conditions(:,1)) = stats.std;
                        
                        SLerrorbar(1:1:360,y,'yError',...
                            ystd,['color=[',num2str(colors{j}),']'],...
                            'Symbol=-','Linestyle','none');
                        bar(1:1:360,y,'Facecolor',colors{j},'edge','none')
                        
                        %motion directions
                        %to the left then to the right
                        dirleftandright(j) = unique(dir.diStimStrengthPrior(C.dispStimStrengthPrior));
                        if ~isnan(dirleftandright(j))
                            plot([dirleftandright(j) dirleftandright(j)],[0 0.5],'-','color',colormk(j,:),'linewidth',3);
                        end
                        
                        %prior mean
                        plot([prior.m prior.m],[0 0.5],'-','color',[0.3 0.6 0.9],'linewidth',3);
                        axis on
                    end
                end
                
                %put on top
                %uistack(h1,'top')
                %uistack(h2,'top')
                
                %label
                set(gca,'xtick',[],'xticklabel',[])
                
                %title the
                if i==1
                    title([num2str(StimStrengthi*100),'% Stim strength and prior std=',num2str(pstdi)],'fontsize',11)
                end
                
                %label last axis
                if i==dist2Pr.uniqCount
                    ylabel(subp(i),'Reaction times','fontsize',11)
                    xlabel({'Estimated directions','(degrees)'},'fontsize',11)
                    set(gca,'xtick',[65 225 385],'xticklabel',[65 225 25])
                end
                xlim([0 360])
                ylim([0 max(pmax(:))])
                %drawPublishAxis
            end
        end
        
        %draw data sorted by signed distance of displayed direction to prior mean
        if strcmp(type,'signDistance')
            
            %plot color
            colors=[1 0 0];
            
            %set of experimental priors
            Thepriors=unique(pstd);
            
            %set of stim strength
            TheStimStrengths=unique(StimStrength);
            
            %motion direction as signed linear distance to prior mean
            dir.dilinDisttoPrior=dir.di-prior.m;
            
            %set of motion directions as signed linear distance to prior mean
            ThedirlinDisttoPrior=unique(dir.dilinDisttoPrior);
            
            %look at how data distribution changes as prior strength increases
            %get number of axes
            axesPos = 1 : 1:numel(ThedirlinDisttoPrior)*numel(Thepriors);
            dirpos = [];
            
            %calculate position of plot for each motion direction
            for i=1:numel(Thepriors)
                dirpos = [dirpos; i:numel(Thepriors):numel(axesPos)];
            end
            
            %look at this stim strength
            for j=1:numel(TheStimStrengths)
                thisStimStrength=TheStimStrengths(j);
                
                %open figure
                scrsz=get(0,'ScreenSize');
                width=7;
                figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
                
                %look at this prior
                for i = 1 : numel(Thepriors)
                    thisprior=Thepriors(i);
                    
                    %Look at this motion direction
                    for k = 1 : numel(ThedirlinDisttoPrior)
                        
                        %get this motion direction
                        %in term of distance to the prior mean
                        thisDirlinDisttoPrior = ThedirlinDisttoPrior(k);
                        
                        %subplots
                        subplot(numel(ThedirlinDisttoPrior),numel(Thepriors),dirpos(i,k));
                        
                        %get data for this condition
                        %dirtoHist=dir.es(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        RTtoHist = databank.reactionTime(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        estimates = databank.estimatesDeg(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        estimates(estimates==0) = 360;
                        
                        %if this condition exists
                        if ~isempty(RTtoHist) && sum(isnan(RTtoHist))~=numel(RTtoHist)
                            
                            %stats
                            stats = SLmakeStat(RTtoHist,estimates);
                            pmax(j) = max(stats.mean);
                            
                            %show sample size
                            RTCount = sum(~isnan(RTtoHist));
                            fprintf('%s %12g \n', '(analyses) number of trial is',RTCount)
                            
                            %data to plot
                            hold all
                            y = nan(360,1);
                            y(stats.conditions(:,1)) = stats.mean;
                            
                            ystd = nan(360,1);
                            ystd(stats.conditions(:,1)) = stats.std;
                            
                            
                            %----
                            %draw
                            %-----
                            SLerrorbar(1:1:360,y,'yError',...
                                ystd,['color=[',num2str(colors),']'],...
                                'Symbol=-','Linestyle','none');
                            bar(1:1:360,y,'Facecolor',colors,'edge','none')
                            
                            %graphics
                            xlim([0 360])
                            ylim([0 max(y)])
                            axis tight
                            set(gca,'ytick',[0 fix(max(y)*100)/100],'ytickLabel',[0 fix(max(y)*100)/100])
                            
                            %indicate motion direction on the plot
                            %get motion direction
                            thisDir = thisDirlinDisttoPrior+prior.m;
                            plot([thisDir thisDir],[0 0.5*max(y)],'-','color',[0 0 0],...
                                'linewidth',3)
                            
                            %prior mean
                            plot([prior.m prior.m],[0 0.5*max(y)],'-','color',...
                                [0.3 0.6 0.9],'linewidth',3);
                            
                            %remove axis for clarity
                            axis on
                            %set(gca,'ytick',0,...
                            %    'xticklabel',0)
                        else
                            axis off
                        end
                        box off
                        
                        %graphics
                        %--------
                        %title top axes
                        if k==1
                            title(['(',num2str(thisStimStrength*100),'% stim strength, std of prior=',num2str(thisprior),')'],'fontsize',11)
                        end
                        %no x-tick
                        set(gca,'xtick',[],'xtickLabel',[],'fontsize',10)
                        
                        %x-axis only for the bottom axis
                        if k==numel(ThedirlinDisttoPrior) && i==1
                            ylabel('Probability','fontsize',11)
                        end
                        
                        %ylabel for one of the axis only
                        if k==numel(ThedirlinDisttoPrior)
                            axis on
                            xlabel('Estimated directions (degrees)','fontsize',11)
                            xlim([0 360])
                            set(gca,'xtick',[65 225 385],'xtickLabel',[65 225 25])
                        end
                    end
                end
            end
        end
    end
end

%%%Reaction times near prior vs near motion direction
function drawRTPriorVsMotDirGrp(databank,type,varargin)

%experimental factors
%motion directions, prior, StimStrength
x = 5 : 10 : 360;
dir.di = cell2mat(databank.data(:,strcmp(databank.nm,'sample_dir')==1));
pstd = cell2mat(databank.data(:,strcmp(databank.nm,'Pstd')==1));
StimStrength = cell2mat(databank.data(:,strcmp(databank.nm,'StimStrength')==1));
prior.m = unique(cell2mat(databank.data(:,strcmp(databank.nm,'priormean')==1)));

%(case von Mises prior)
%---------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
        
        %draw data sorted by signed distance of displayed direction to prior mean
        if strcmp(type,'signDistance')
            
            %plot color
            colors = [.5 .5 .5];
            
            %set of experimental priors
            Thepriors = unique(pstd);
            
            %set of stim strength
            TheStimStrengths = unique(StimStrength);
            
            %motion direction as signed linear distance to prior mean
            dir.dilinDisttoPrior = dir.di - prior.m;
            
            %set of motion directions as signed linear distance to prior mean
            ThedirlinDisttoPrior = unique(dir.dilinDisttoPrior);
            
            %look at how data distribution changes as prior strength increases
            %get number of axes
            axesPos = 1 : 1 : numel(ThedirlinDisttoPrior)*numel(Thepriors);
            dirpos = [];
            
            %calculate position of plot for each motion direction
            for i = 1 : numel(Thepriors)
                dirpos = [dirpos; i:numel(Thepriors):numel(axesPos)];
            end
            
            %!!!!!!!!!!!!!!test code!!!!!!!!!!!!!!!
            %uncomment to use
            %subject randomly report 225 or displayed direction
            %we set RT higher for for prior estimate do we see
            %that in the graph? Yes as expected.
            %databank.estimatesDeg = (rand(numel(databank.estimatesDeg),1)>0.5)*225;
            %databank.estimatesDeg(databank.estimatesDeg==0) = dir.di(databank.estimatesDeg==0);
            %databank.reactionTime(databank.estimatesDeg==225)=100;
            
            %look at this stim strength
            for j = 1 : numel(TheStimStrengths)
                thisStimStrength = TheStimStrengths(j);
                
                %open figure
                scrsz=get(0,'ScreenSize');
                width=7;
                figure('color','w','Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
                
                %look at this prior
                for i = 1 : numel(Thepriors)
                    thisprior = Thepriors(i);
                    
                    %Look at this motion direction
                    for k = 1 : numel(ThedirlinDisttoPrior)
                        
                        %get this motion direction
                        %in term of distance to the prior mean
                        thisDirlinDisttoPrior = ThedirlinDisttoPrior(k);
                        
                        %subplots
                        subplot(numel(ThedirlinDisttoPrior),numel(Thepriors),dirpos(i,k));
                        
                        %get data for this condition
                        %dirtoHist=dir.es(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        RTtoHist = databank.reactionTime(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        estimates = databank.estimatesDeg(pstd==thisprior & StimStrength==thisStimStrength & dir.dilinDisttoPrior==thisDirlinDisttoPrior);
                        estimates(estimates==0) = 360;
                        
                        %if this condition exists and displayed direction
                        %is not at the prior
                        if ~isempty(RTtoHist) && sum(isnan(RTtoHist))~=numel(RTtoHist) && thisDirlinDisttoPrior~=0
                            
                            %group estimates by distance to prior and motion
                            %direction
                            distEs2PrThisC = abs(SLvectors2signedAngle(estimates,prior.m));
                            distEs2DiThisC = abs(SLvectors2signedAngle(estimates,thisDirlinDisttoPrior + prior.m));
                            
                            %RT clustered for estimates near prior (2) vs near direction (1)
                            ClustIdx = [distEs2PrThisC < distEs2DiThisC]+1;
                            RTstats = SLmakeStat(RTtoHist,ClustIdx);
                            pmax(j) = max(RTstats.mean);
                            
                            %sample size
                            RTcount = nan(1,2);
                            RTcount(RTstats.conditions(:,1)) = RTstats.count;
                            fprintf('%s \n',['(analyses) number of trials are:' num2str(RTcount)])
                            
                            %----
                            %draw
                            %-----
                            hold all
                            bar(RTstats.conditions(:,1),RTstats.mean,'Facecolor',colors,'edge','none')
                            
                            SLerrorbar(RTstats.conditions(:,1),RTstats.mean,'yError',...
                                RTstats.std,['color=[',num2str(colors),']'],...
                                'Symbol=-','Linestyle','none');
                            
                            
                            %graphics
                            xlim([0 3])
                            ylim([0 max(RTstats.mean)])
                            set(gca,'ytick',[0 fix(max(RTstats.mean)*100)/100],...
                                'ytickLabel',[0 fix(max(RTstats.mean)*100)/100])
                            
                            %indicate motion direction on the plot
                            %get motion direction
                            thisDir = thisDirlinDisttoPrior + prior.m;
                            plot([thisDir thisDir],[0 0.5*max(RTstats.mean)],'-','color',[0 0 0],...
                                'linewidth',3)
                            
                            %prior mean
                            plot([prior.m prior.m],[0 0.5*max(RTstats.mean)],'-','color',...
                                [0.3 0.6 0.9],'linewidth',3);
                            
                            %remove axis for clarity
                            axis on
                        else
                            axis off
                        end
                        box off
                        
                        %graphics
                        %--------
                        %title top axes
                        if k==1
                            title(['(',num2str(thisStimStrength*100),'% StimStrengtherence, std of prior=',num2str(thisprior),')'],'fontsize',11)
                        end
                        %no x-tick
                        set(gca,'xtick',[],'xtickLabel',[],'fontsize',10)
                        
                        %x-axis only for the bottom axis
                        if k==numel(ThedirlinDisttoPrior) && i==1
                            ylabel('Probability','fontsize',11)
                        end
                        
                        %ylabel for one of the axis only
                        if k==numel(ThedirlinDisttoPrior)
                            axis on
                            xlabel({'Estimates clustered into' ,'"near prior" or "near motion direction" (degrees)'},'fontsize',11)
                            xlim([0 3])
                            %                             set(gca,'xtick',25:40:355,'xtickLabel',25:40:355)
                        end
                    end
                end
            end
        end
    end
end

%%raw data
function drawRaw(databank)

%data
data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'estimatedFeature'))==1)));

%make sure there is a unique value 360 for 0 and 360
data(data==0)=360;

%get the factors
d=cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));
coh=cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
pstd=cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));

%factors 1,2,3
F.f1.i=d;
F.f1.nm='d';
F.f1.L=unique(F.f1.i);
F.f1.L=sort(F.f1.L,'ascend');
F.f1.n=numel(F.f1.L);

F.f2.i=StimStrength;
F.f2.nm='StimStrength';
F.f2.L=unique(F.f2.i);
F.f2.L=sort(F.f2.L,'descend');
F.f2.n=numel(F.f2.L);

F.f3.i=pstd;
F.f3.nm='pstd';
F.f3.L=unique(F.f3.i);
F.f3.L=sort(F.f3.L,'descend');
F.f3.n=numel(F.f3.L);

%positions main
for i=1:F.f1.n
    F.f1.pos(i)={find(F.f1.i==F.f1.L(i))};
end
for i=1:F.f2.n
    F.f2.pos(i)={find(F.f2.i==F.f2.L(i))};
end
for i=1:F.f3.n
    F.f3.pos(i)={find(F.f3.i==F.f3.L(i))};
end

%positions inter
for k=1:F.f1.n
    for j=1:F.f2.n
        for i=1:F.f3.n
            F.inter.pos(k,i,j)=...
                {intersect( ...
                intersect(F.f1.pos{k},F.f2.pos{j}),...
                F.f3.pos{i})};
        end
    end
end

%Graphics
fg2=figure('color','w');
set(fg2, 'units', 'centimeters', 'pos', [0 0 30 10])
F.f2.color0=[0.3 0 0;...
    1 0.2 0;...
    1 0.6 0;...
    0.75 0.75 0];
one=ones(size(F.f2.color0,1),3);

%marker's edges, color & size
egdec=F.f2.color0+1*(one-F.f2.color0);
egdecrw=F.f2.color0+0.3*(one-F.f2.color0);
zer=zeros(size(F.f2.color0,1),3);
colPred=F.f2.color0+0.3*(zer-F.f2.color0);
F.f2.color=F.f2.color0+0.5*(one-F.f2.color0);
markersizeS=40;
ft=14;

for j=1:F.f2.n
    hs(1)=subplot(1,3,j);
    axis square
    for k=1:F.f1.n
        for i=1:F.f3.n
            rawData{k,i,j}=data(F.inter.pos{k,i,j},:);
            
            %data
            hold all
            subplot(hs(1))
            scatter(hs(1),d(F.inter.pos{k,i,j}),rawData{k,i,j}',...
                'MarkerEdgeColor',egdecrw(i,:),...
                'MarkerFaceColor',F.f2.color0(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            
            %ideal estimation
            subplot(hs(1))
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            
            %indicate the prior mean if the experiment is with Gaussian
            %prior
            plot(F.f1.L,225*ones(F.f1.n,1),'b.','Markersize',4)
            
            %graphics
            %--------
            xlim([0 360])
            ylim([0 360])
            ylabel(hs(1),{'Estimated directions', '(degrees)'},'fontsize',ft)
            if j==2
                xlabel(hs(1),{'Displayed directions', '(degrees)'},'fontsize',ft)
            end
        end
    end
    axis tight
end

%%Test "sequential effect" hypothesis
function dataPastEffect01(databank,FA,varargin)
%The data could be explained by two mechanisms (long or short-term memory)
%
%   Prior effect (trial history)
%   ----------------------------
%   Subjects learnt the mean of the experimental prior and
%   use it to guess the noisy motion direction in the current trial.
%   Bayesian and Competition mechanisms are both consistent with this
%   mechanism as those mechanism represent and use the prior
%   (mean and strength).
%
%   Recency effect (just previous trial)
%   ------------------------------------
%   Subjects are simply biased by the past trial's motion
%   direction and do not learn the prior.
%
%So, are data attracted toward the mean the experimental prior (prior
%effect) or toward the past trial's motion direction (recency effect)?

fprintf('%s \n','(analyses) Now testing for "sequential effect"...')

if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
        
        %Test the code (all tests worked)
        %----------------------------------------------------------------------
        %%simulate when subjects estimate the current direction: arrow should point
        %current direction if code correct
        %dir.estiCoor=polar2cartesian(FA.g3.thisT,1);
        
        %%simulate when subjects estimate the current direction: arrow should point
        %current direction if code correct
        %dir.estiCoor=polar2cartesian(FA.g3.thisT,1);
        
        %%simulate when subjects estimate the previous direction: arrow should point
        %previous direction if code correct
        %dir.estiCoor=polar2cartesian([NaN;FA.g3.thisT(1:end-1)],1);
        
        %%simulate when subjects estimate the prior mean: arrow should point prior
        %mean if code correct
        %dir.estiCoor=polar2cartesian(225*ones(size(databank.data,1),1),1);
        %----------------------------------------------------------------------
        
        %test the code (works fine)
        %--------------------------------
        %Here simulated estimate should point toward previous direction
        %simulate choosing previous displayed
        %         dir.estiCoor=repmat(polar2cartesian(135,1),length(FA.g3.thisT),1);
        %
        %         %displayed alternates between 135 and 180
        %         dir.dispCoor=repmat(polar2cartesian([135;180],[1 1]),length(FA.g3.thisT),1);
        %
        %         %Get the coordinates of the prior mean
        %         dir.priormeanCoor=polar2cartesian(225,1);
        
        %(Case von Mises prior experiment)
        %--------------------------------
        %Get the coordinates of the estimated directions (actual code)
        dir.estiCoor=cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
        
        %Get the coordinates of the displayed directions
        r=1;
        dir.dispCoor=polar2cartesian(FA.g3.thisT,r);
        
        %Get the coordinates of the prior mean
        dir.priormeanCoor=polar2cartesian(225,r);
        
        %Get the directions that were displayed between the mean of the prior and
        %the past trial's directions.
        %e.g., 225? (prior mean) < 180? (current) < 135? (previous)
        %Calculate the average direction displayed on current trials, on past
        %trial and look at where average estimate in current trials lean
        %towards: prior or past trial?
        %Look at each factor 1
        meanDispDir=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        meanPrevDir=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        meanEstimate=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        
        %calculate and draw
        for j=1:FA.g1.lvlsnb
            
            %find this factors' trials
            trialsThisFA1=FA.g1.thisT==FA.g1.lvlsnm(j);
            
            %look at factor 2
            for k=1:FA.g2.lvlsnb
                
                %find this factor's trials
                trialsThisFA2=FA.g2.thisT==FA.g2.lvlsnm(k);
                
                %get positions of current trials with this condition
                posCurrentThisCondition=find(trialsThisFA1 & trialsThisFA2);
                
                %get positions of the trials that preceded
                posPrevious=posCurrentThisCondition-1;
                
                %remove first trials because there is not preceding trial
                posPrevious(posCurrentThisCondition==1)=[];
                posCurrentThisCondition(posCurrentThisCondition==1)=[];
                
                %get current and previous direction
                dir.dispCoorCurThisCondition=dir.dispCoor(posCurrentThisCondition,:);
                dir.dispCoorPrev=dir.dispCoor(posPrevious,:);
                
                %get the current directions displayed in our condition
                for i=1:size(dir.dispCoorCurThisCondition,1)
                    
                    %When current and previous directions were displayed
                    %counterclockwise to prior's mean (i.e.,0< angle(prior,current)<180)
                    angle.PriortoCur(i)=SLvectors2signedAngle(dir.priormeanCoor,...
                        dir.dispCoorCurThisCondition(i,:));
                    
                    %angle between current and previous direction
                    angle.CurtoPrev(i)=SLvectors2signedAngle(dir.dispCoorCurThisCondition(i,:),...
                        dir.dispCoorPrev(i,:));
                    
                    %get current trials when the sequence is: past < current < prior
                    %case current counterclockwise to prior mean
                    %case previous counterclockwise to current for this condition
                    if  0<angle.PriortoCur(i) && angle.PriortoCur(i)<180 && 0<angle.CurtoPrev(i) && angle.CurtoPrev(i)<180
                        
                        %store the position of those trials
                        trialCounterClock{j,k}(i)=posCurrentThisCondition(i);
                    else
                        trialCounterClock{j,k}(i)=NaN;
                    end
                end
                
                %now get rid of NaN
                trialCounterClock{j,k}(isnan(trialCounterClock{j,k}))=[];
                
                %averages current, previous directions & estimates
                %-------------------------------------------------
                %displayed direction in current trial
                meanDispDirInfo{j,k}=vectorStat(dir.dispCoor(trialCounterClock{j,k},:));
                meanDispDir{j,k}=meanDispDirInfo{j,k}.coord.mean;
                
                %estimated direction in current trial
                meanEstimateInfo{j,k}=vectorStat(dir.estiCoor(trialCounterClock{j,k},:));
                meanEstimate{j,k}=meanEstimateInfo{j,k}.coord.mean;
                
                %displayed direction in previous trial
                meanPrevDirInfo{j,k}=vectorStat(dir.dispCoor(trialCounterClock{j,k}-1,:));
                meanPrevDir{j,k}=meanPrevDirInfo{j,k}.coord.mean;
            end
        end
        
        %draw all conditions
        meanPriorCoord=polar2cartesian(225,10);
        count=0;
        
        %graphics
        fig=figure('color','w');
        ax1=[];
        FA1nmlabel=[];
        FA2nmlabel=[];
        
        %look at each condition
        for k=1:FA.g2.lvlsnb
            
            for j=1:FA.g1.lvlsnb
                
                %draw Estimates density
                %----------------------
                %get all estimates of current trials directions in this condition
                CurEstimateThisCondition=meanEstimateInfo{j,k}.deg.all;
                
                %get number of estimates
                numCurEstimateThisCondition=length(CurEstimateThisCondition);
                
                %count occurrence of each estimated direction
                x=5:10:355;
                probaCurEstimateThisCondition=histc(CurEstimateThisCondition,x)/numCurEstimateThisCondition;
                
                %initialize axes and axes' legends
                count=count+1;
                ax1=[ax1 subplot(numel(FA.g2.lvlsnm),numel(FA.g1.lvlsnm),count)];
                FA1nmlabel=[FA1nmlabel FA.g1.lvlsnm(j)];
                FA2nmlabel=[FA2nmlabel FA.g2.lvlsnm(k)];
                
                %get the maximum of the estimate density to scale mean estimate
                %vector
                maxData=max(probaCurEstimateThisCondition);
                
                %increase the polar radius to get a bit of space for visibility
                radius=1.1*maxData;
                %SLpolar(de2r(x,[]),probaCurEstimateThisCondition',[0.7 0 0],'xSpokes',...
                %    'SpokesTickSteps',8,'area',...
                %    'facecolor',[1 0.5 0.5],...
                %    'edgecolor',[1 0.5 0.5],...
                %    'linestyle','-')
                
                SLpolar(de2r(x,[]),probaCurEstimateThisCondition',[0.7 0 0],...
                    'edgecolor',[1 0.5 0.5],...
                    'linestyle','-')
                
                
                %draw the mean estimate vector
                %-----------------------------
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanEstimate{j,k}(1),meanEstimate{j,k}(2));
                
                %this scales the arrow to the size of the radius
                putOncircle=Vectnorm(j,k)/maxData;
                
                %draw mean estimated direction
                hold on
                arrow([0 0],[meanEstimate{j,k}(1) meanEstimate{j,k}(2)]/putOncircle,...
                    'Length',2,'BaseAngle',[],'Width',1,...
                    'facecolor',[0.5 0 0],...
                    'edgecolor',[0.5 0 0],...
                    'linestyle','-')
                
                %                 %draw previous direction density
                %                 %------------------------------
                %                 %get all estimates of current trials directions in this condition
                %                 prevDirThisCondition=meanPrevDirInfo{j,k}.deg.all;
                %
                %                 %get number of estimates
                %                 numprevDirThisCondition=length(prevDirThisCondition);
                %
                %                 %count occurrence of each estimated direction
                %                 x=5:10:355;
                %                 probaprevDirThisCondition=histc(prevDirThisCondition,x)/numprevDirThisCondition;
                %
                %                 %get the maximum of the previous displayed direction
                %                 %density to scale mean estimate vector
                %                 maxData=max(probaprevDirThisCondition);
                %
                %                 %increase the polar radius to get a bit of space for visibility
                %                 radius=1.1*maxData;
                %                 SLpolar(de2r(x,[]),probaprevDirThisCondition',[0.7 0 0],'xSpokes',...
                %                     'SpokesTickSteps',4,...
                %                     'edgecolor',[.5 .5 .5],...
                %                     'linestyle','--')
                
                %draw current direction
                %----------------------
                %indicate the mean of the current direction
                %get all current trials displayed directions in this condition
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanDispDir{j,k}(1),meanDispDir{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanDispDir{j,k}(1)/putOncircle,meanDispDir{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[.2 .2 .2])
                
                %draw previous direction
                %-----------------------
                %indicate the mean of the current direction
                %get all current trials displayed directions in this condition
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanPrevDir{j,k}(1),meanPrevDir{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanPrevDir{j,k}(1)/putOncircle,meanPrevDir{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[.7 .7 .7])
                
                
                %draw prior mean
                %---------------
                %indicate the mean of the prior
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanPriorCoord(1),meanPriorCoord(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanPriorCoord(1)/putOncircle,meanPriorCoord(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[0 0.7 1])
            end
        end
        
        %legend the axes
        for i=1:numel(ax1)
            title(ax1(i),[num2str(100*FA1nmlabel(i)),'% ',FA.g1.nm,' - ',...
                num2str(FA2nmlabel(i)),' degrees ',FA.g2.nm],'fontsize',8,...
                'fontweight','Bold')
        end
        
        %remove dead spaces
        %SLremoveDeadSpace
        
    end
end


%(Case bimodal prior experiment)
%-------------------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')%
        
        %Get the coordinates of the estimated directions (actual code)
        dir.estiCoor=cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
        
        %Test the code (all tests worked)
        %----------------------------------------------------------------------
        %%simulate when subjects estimate the current direction: arrow should point
        %current direction if code correct
        %dir.estiCoor=polar2cartesian(FA.g3.thisT,1);
        
        %%simulate when subjects estimate the previous direction: arrow should point
        %previous direction if code correct
        %dir.estiCoor=polar2cartesian([NaN;FA.g3.thisT(1:end-1)],1);
        %----------------------------------------------------------------------
        
        %set the radius of the circle but the true size is 2.5 degrees
        %(just for visualization purpose)
        r=1;
        
        %Get the coordinates of the displayed directions
        dir.dispCoor=polar2cartesian(FA.g3.thisT,r);
        
        %counterclock mode at each trial
        %when the distance between the first and second mode is <0 the
        %first mode is more counterclockwise than the second and when it is
        %>0, the second mode is more clockwise than the first. Now get the
        %coordinates of the counterclock mode (left)
        distModes=SLvectors2signedAngle(FA.g2.thisT(:,1),FA.g2.thisT(:,2));
        cclockMode(distModes<0)=FA.g2.thisT(distModes<0,1);
        cclockMode(distModes>0)=FA.g2.thisT(distModes>0,2);
        dir.priorcclockModeCoor=polar2cartesian(cclockMode,r);
        
        %Now get the coordinates of the clock mode (right)
        clockMode(distModes<0)=FA.g2.thisT(distModes<0,2);
        clockMode(distModes>0)=FA.g2.thisT(distModes>0,1);
        dir.priorclockModeCoor=polar2cartesian(clockMode,r);
        
        %Get the directions that were displayed between the mean of the prior and
        %the past trial's directions.
        %e.g., 225? (prior mean) < 180? (current) < 135? (previous)
        %Calculate teh average direction displayed on current trials, on past
        %trial and look at where average estimate in current trials lean
        %towards: prior or past trial?
        %Look at each factor 1
        meanDispDir=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        meanPrevDir=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        meanEstimate=num2cell(nan(FA.g1.lvlsnb,FA.g2.lvlsnb));
        
        %case factor 1 is the prior modes
        if strcmp(FA.g1.nm,'priormodes')
            FA.g1.thisTSeqAna=FA.g1.distance_mode;
            %case not prior mode
        else
            FA.g1.thisTSeqAna=FA.g1.thisT;
        end
        
        %case factor 2 is the prior modes
        if strcmp(FA.g2.nm,'priormodes')
            FA.g2.thisTSeqAna=FA.g2.distance_mode;
            %case not prior mode
        else
            FA.g2.thisTSeqAna=FA.g2.thisT;
        end
        
        %case factor 3 is the prior modes
        if strcmp(FA.g3.nm,'priormodes')
            FA.g3.thisTSeqAna=FA.g3.distance_mode;
            %case not prior mode
        else
            FA.g3.thisTSeqAna=FA.g3.thisT;
        end
        
        %calculate and draw
        for j=1:FA.g1.lvlsnb
            
            %find this factors' trials
            trialsThisFA1=FA.g1.thisTSeqAna==FA.g1.lvlsnm(j);
            
            %look at factor 2
            for k=1:FA.g2.lvlsnb
                
                %find this factor's trials
                trialsThisFA2=FA.g2.thisTSeqAna==FA.g2.lvlsnm(k);
                
                %get positions of current trials with this condition
                posCurrentThisCondition=find(trialsThisFA1 & trialsThisFA2);
                
                %get positions of the trials that preceded
                posPrevious=posCurrentThisCondition-1;
                
                %remove first trials because there is not a preceding trial
                posPrevious(posCurrentThisCondition==1)=[];
                posCurrentThisCondition(posCurrentThisCondition==1)=[];
                
                %get current trial's motion direction
                dir.dispCoorCurThisCondition=dir.dispCoor(posCurrentThisCondition,:);
                
                %get counterclock mode in the current trial
                dir.priorcclockModeCoorThisCondition=dir.priorcclockModeCoor(posCurrentThisCondition,:);
                
                %get the previous motion direction
                dir.dispCoorPrev=dir.dispCoor(posPrevious,:);
                
                
                
                %get the current directions
                for i=1:size(dir.dispCoorCurThisCondition,1)
                    
                    %When current and previous directions were displayed
                    %counterclockwise to the most counterclockwise mode of
                    %the bimodal prior (i.e.,0< angle(left prior mode,current)<180)
                    angle.PriortoCur(i)=SLvectors2signedAngle(dir.priorcclockModeCoorThisCondition(i,:),...
                        dir.dispCoorCurThisCondition(i,:));
                    
                    %angle between current and previous direction
                    angle.CurtoPrev(i)=SLvectors2signedAngle(dir.dispCoorCurThisCondition(i,:),...
                        dir.dispCoorPrev(i,:));
                    
                    %get current trials when the sequence is: past < current < prior
                    %case current counterclockwise to prior mean
                    %case previous counterclockwise to current for this condition
                    if  0<angle.PriortoCur(i) && angle.PriortoCur(i)<180 && 0<angle.CurtoPrev(i) && angle.CurtoPrev(i)<180
                        
                        %store the position of those trials
                        trialCounterClock{j,k}(i)=posCurrentThisCondition(i);
                    else
                        trialCounterClock{j,k}(i)=NaN;
                    end
                end
                
                %now get rid of NaN
                trialCounterClock{j,k}(isnan(trialCounterClock{j,k}))=[];
                
                %averages current, previous directions & estimates
                %-------------------------------------------------
                %displayed direction in current trial
                meanDispDirInfo{j,k}=vectorStat(dir.dispCoor(trialCounterClock{j,k},:));
                meanDispDir{j,k}=meanDispDirInfo{j,k}.coord.mean;
                
                %estimated direction in current trial
                meanEstimateInfo{j,k}=vectorStat(dir.estiCoor(trialCounterClock{j,k},:));
                meanEstimate{j,k}=meanEstimateInfo{j,k}.coord.mean;
                
                %counterclock prior mode in current trial
                meancclockPmodeInfo{j,k}=vectorStat(dir.priorcclockModeCoor(trialCounterClock{j,k},:));
                meancclockPmode{j,k}=meancclockPmodeInfo{j,k}.coord.mean;
                
                %clock prior mode in current trial
                meanclockPmodeInfo{j,k}=vectorStat(dir.priorclockModeCoor(trialCounterClock{j,k},:));
                meanclockPmode{j,k}=meanclockPmodeInfo{j,k}.coord.mean;
                
                %displayed direction in previous trial
                meanPrevDirInfo{j,k}=vectorStat(dir.dispCoor(trialCounterClock{j,k}-1,:));
                meanPrevDir{j,k}=meanPrevDirInfo{j,k}.coord.mean;
            end
        end
        
        %draw all conditions
        count=0;
        
        %graphics
        figure('color','w');
        ax1=[];
        FA1nmlabel=[];
        FA2nmlabel=[];
        
        %look at each condition
        for k=1:FA.g2.lvlsnb
            
            for j=1:FA.g1.lvlsnb
                
                %draw Estimates density
                %----------------------
                %get all estimates of current trials directions in this condition
                CurEstimateThisCondition=meanEstimateInfo{j,k}.deg.all;
                
                %get number of estimates
                numCurEstimateThisCondition=length(CurEstimateThisCondition);
                
                %count occurrence of each estimated direction
                x=5:10:355;
                probaCurEstimateThisCondition=histc(CurEstimateThisCondition,x)/numCurEstimateThisCondition;
                
                %initialize axes and axes' legends
                count=count+1;
                ax1=[ax1 subplot(numel(FA.g2.lvlsnm),numel(FA.g1.lvlsnm),count)];
                FA1nmlabel=[FA1nmlabel FA.g1.lvlsnm(j)];
                FA2nmlabel=[FA2nmlabel FA.g2.lvlsnm(k)];
                
                %get the maximum of the estimate density to scale mean estimate
                %vector
                maxData=max(probaCurEstimateThisCondition);
                
                %increase the polar radius to get a bit of space for visibility
                radius=1.1*maxData;
                
                %if we found data when motion direction is in between counterclock
                %prior mode and previous trial
                if ~isnan(CurEstimateThisCondition)
                    
                    %plot polar
                    SLpolar(de2r(x,[]),probaCurEstimateThisCondition,[0.7 0 0],'xSpokes',...
                        'SpokesTickSteps',4,'area',...
                        'facecolor',[1 0.5 0.5],...
                        'edgecolor',[1 0.5 0.5])
                    
                    %draw the mean estimate vector
                    %get lenth of vector
                    Vectnorm(j,k)=SLcalculateVectNorm(meanEstimate{j,k}(1),meanEstimate{j,k}(2));
                    
                    %this scales the arrow to the size of the radius
                    putOncircle=Vectnorm(j,k)/maxData;
                    
                else
                    %plot blank polar and null variable
                    SLpolar(de2r(x,[]),0.1*ones(numel(de2r(x,[])),1),[1 1 1],'xSpokes',...
                        'SpokesTickSteps',4,'area',...
                        'facecolor',[1 1 1],...
                        'edgecolor',[1 1 1])
                    meanEstimate{j,k}=[0 0];
                    putOncircle=1;
                end
                
                %draw mean estimated direction
                hold on
                arrow([0 0],[meanEstimate{j,k}(1) meanEstimate{j,k}(2)]/putOncircle,...
                    'Length',2,'BaseAngle',[],'Width',1,...
                    'facecolor',[0.5 0 0],...
                    'edgecolor',[0.5 0 0])
                
                %draw current direction
                %----------------------
                %indicate the mean of the current direction
                %get all current trials displayed directions in this condition
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanDispDir{j,k}(1),meanDispDir{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanDispDir{j,k}(1)/putOncircle,meanDispDir{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[.2 .2 .2])
                
                %draw previous direction
                %-----------------------
                %indicate the mean of the current direction
                %get all current trials displayed directions in this condition
                %get lenth of vector
                Vectnorm(j,k)=SLcalculateVectNorm(meanPrevDir{j,k}(1),meanPrevDir{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnorm(j,k)/radius;
                hold on
                plot(meanPrevDir{j,k}(1)/putOncircle,meanPrevDir{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[.7 .7 .7])
                
                
                %draw average prior modes
                %------------------------
                %indicate the counterclock mode of the prior
                %get lenth of vector
                Vectnormcclock(j,k)=SLcalculateVectNorm(meancclockPmode{j,k}(1),meancclockPmode{j,k}(2));
                
                %scale its length to the size of the radius
                putOncircle=Vectnormcclock(j,k)/radius;
                hold on
                plot(meancclockPmode{j,k}(1)/putOncircle,meancclockPmode{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[0 0.7 1])
                
                %indicate the clock mode of the prior
                Vectnormclock(j,k)=SLcalculateVectNorm(meanclockPmode{j,k}(1),meanclockPmode{j,k}(2));
                putOncircle=Vectnormclock(j,k)/radius;
                hold on
                plot(meanclockPmode{j,k}(1)/putOncircle,meanclockPmode{j,k}(2)/putOncircle,'.',...
                    'MarkerSize',20,...
                    'color',[0 0.7 1])
            end
        end
        
        %legend the axes
        for i=1:numel(ax1)
            title(ax1(i),[num2str(100*FA1nmlabel(i)),'% ',FA.g1.nm,' - ',...
                num2str(FA2nmlabel(i)),' degrees (mode1-mode2) ',FA.g2.nm],'fontsize',8,...
                'fontweight','Bold')
        end
        
        %remove dead spaces
        %SLremoveDeadSpace(0.05)
    end
end

%%Other way to test "sequential effect" hypothesis
function dataPastEffect02(databank,FA,varargin)

fprintf('%s \n','(dataPastEffect02) Now testing for "sequential effect"...')

%(Case von Mises prior experiment)
%-------------------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')
        
        %Test the code (all tests worked)
        %----------------------------------------------------------------------
        %%simulate when subjects estimate the current direction: arrow should point
        %current direction if code correct
        %dir.estiCoor=polar2cartesian(FA.g3.thisT,1);
        
        %%simulate when subjects estimate the current direction: arrow should point
        %current direction if code correct
        %dir.estiCoor=polar2cartesian(FA.g3.thisT,1);
        
        %%simulate when subjects estimate the previous direction: arrow should point
        %previous direction if code correct
        %dir.estiCoor=polar2cartesian([NaN;FA.g3.thisT(1:end-1)],1);
        
        %%simulate when subjects estimate the prior mean: arrow should point prior
        %mean if code correct
        %dir.estiCoor=polar2cartesian(225*ones(size(databank.data,1),1),1);
        %----------------------------------------------------------------------
        
        %test the code (works fine)
        %--------------------------------
        %Here simulated estimate should point toward previous direction
        %simulate choosing previous displayed
        %         dir.estiCoor=repmat(polar2cartesian(135,1),length(FA.g3.thisT),1);
        %
        %         %displayed alternates between 135 and 180
        %         dir.dispCoor=repmat(polar2cartesian([135;180],[1 1]),length(FA.g3.thisT),1);
        %
        %         %Get the coordinates of the prior mean
        %         dir.priormeanCoor=polar2cartesian(225,1);
        
        %Get coordinates
        %prior mean, displayed and estimated directions
        dir.estiCoor=cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
        dir.dispCoor=polar2cartesian(FA.g3.thisT,1);
        dir.priormeanCoor=polar2cartesian(225,1);
        
        %Get trial sequence
        [trialClock,trialCClock]=getTrialSequences(dir,FA,varargin{1});
    end
end

%(Case bimodal prior experiment)
%-------------------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
        
        %Get coordinates
        %prior mean, displaye and estimated directions
        dir.estiCoor=cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
        dir.dispCoor=polar2cartesian(FA.g3.thisT,1);
        
        %Get trial sequence
        [trialClock,trialCClock]=getTrialSequences(dir,FA,varargin{1});
    end
end

%To plot the sequential effect we need to sort trials to make
%sure that
%
%   -previous trial is never the prior mean in which case a bias
%   toward the prior mean will be interpreted as a bias toward the
%   previous trial
%
%   -previous trial is never on the same side as the prior mean with
%   respect to current for the same reason as above.
count=0;
figure('color','w')
for j=1:FA.g1.lvlsnb
    for k=1:FA.g2.lvlsnb
        
        count=count+1;
        subplot(FA.g1.lvlsnb,FA.g2.lvlsnb,count)
        axis square
        
        %Previous is clockwise to current then to prior
        %----------------------------------------------
        %error current: estimate - displayed (deg)
        %angle>0 means that the estimate is more
        %clockwise than the displayed direction
        errorCurTprevClock=SLvectors2signedAngle(dir.dispCoor(trialClock{j,k},:),dir.estiCoor(trialClock{j,k},:));
        
        %relative direction of previous trial (deg)
        %angle>0 means that the previous is more
        %clockwise than the displayed direction
        cur2prevTClock=SLvectors2signedAngle(dir.dispCoor(trialClock{j,k},:),dir.dispCoor(trialClock{j,k}-1,:));
        
        %Previous is counterclockwise to current then to prior
        %-----------------------------------------------------
        errorCurTprevCClock=SLvectors2signedAngle(dir.dispCoor(trialCClock{j,k},:),dir.estiCoor(trialCClock{j,k},:));
        cur2prevTCClock=SLvectors2signedAngle(dir.dispCoor(trialCClock{j,k},:),dir.dispCoor(trialCClock{j,k}-1,:));
        
        %stack them together
        errorCurT{j,k}=[errorCurTprevCClock ; errorCurTprevClock];
        cur2prevT{j,k}=[cur2prevTCClock ; cur2prevTClock];
        
        %check if data are missing
        if isempty(errorCurT{j,k})
            errorCurT{j,k}=NaN;
        end
        if isempty(cur2prevT{j,k})
            cur2prevT{j,k}=NaN;
        end
        
        %draw
        hold all
        scatter(cur2prevT{j,k},errorCurT{j,k},65,[.9 0 0],'markerfacecolor',...
            [.8 0 0])
        
        %title
        title(['StimStrength:',num2str(FA.g1.lvlsnm(j)),' prior:',num2str(FA.g2.lvlsnm(k))])
        
        %graphics
        %--------
        %when data are not missing
        if ~isnan(errorCurT{j,k}) & ~isnan(cur2prevT{j,k})
            
            plot([-180 180],[0 0],'k:')
            plot([0 0],[-180 180],'k:')
            
            %limits
            xlim([-180 180])
            ylim([-180 180])
            set(gca,'xtick',-180:45:180,...
                'xticklabel',-180:45:180)
            
            set(gca,'ytick',-180:45:180,...
                'yticklabel',-180:45:180)
        else
            axis off
        end
        
        %label last axe
        if count==FA.g1.lvlsnb*FA.g2.lvlsnb
            xlabel('Relative direction of previous trial (deg)','fontsize',11)
        end
        
        %label first axe
        ylabpos=FA.g2.lvlsnb*[0:1:FA.g1.lvlsnb-1]+1;
        if sum(count==ylabpos)
            ylabel('Error on current trial (deg)')
        end
        %drawPublishAxis
    end
end
%SLremoveDeadSpace

%%get certain trial sequences (sequential effect)
function [trialClock,trialCClock] = getTrialSequences(dir,FA,varargin)

%(case von Mises prior)
%--------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'vonMisesPrior')%
        %calculate and draw
        for j=1:FA.g1.lvlsnb
            
            %find this factors' trials
            trialsThisFA1=FA.g1.thisT==FA.g1.lvlsnm(j);
            
            %look at factor 2
            for k=1:FA.g2.lvlsnb
                
                %find this factor's trials
                trialsThisFA2=FA.g2.thisT==FA.g2.lvlsnm(k);
                
                %get positions of current trials with this condition
                posCurrentThisCondition=find(trialsThisFA1 & trialsThisFA2);
                
                %get positions of the trials that preceded
                posPrevious=posCurrentThisCondition-1;
                
                %remove first trials because there is not preceding trial
                posPrevious(posCurrentThisCondition==1)=[];
                posCurrentThisCondition(posCurrentThisCondition==1)=[];
                
                %get current and previous direction
                dir.dispCoorCurThisCondition=dir.dispCoor(posCurrentThisCondition,:);
                dir.dispCoorPrev=dir.dispCoor(posPrevious,:);
                
                %get the current directions displayed in our condition
                for i=1:size(dir.dispCoorCurThisCondition,1)
                    
                    %When current and previous directions were displayed
                    %counterclockwise to prior's mean (i.e.,0< angle(prior,current)<180)
                    angle.PriortoCur(i)=SLvectors2signedAngle(dir.priormeanCoor,...
                        dir.dispCoorCurThisCondition(i,:));
                    
                    %angle between current and previous direction
                    angle.CurtoPrev(i)=SLvectors2signedAngle(dir.dispCoorCurThisCondition(i,:),...
                        dir.dispCoorPrev(i,:));
                    
                    %Previous is Clockwise to current then to prior
                    %------------------------------------------
                    %get current trials when the sequence is: prior > current <
                    %previous. We set 179 and not 180 to be sure they are not
                    %at opposite direction in which case interpretation is not
                    %possible
                    if  0<angle.PriortoCur(i) && angle.PriortoCur(i)<179 && 0<angle.CurtoPrev(i) && angle.CurtoPrev(i)<179
                        
                        %store the position of those trials
                        trialClock{j,k}(i)=posCurrentThisCondition(i);
                    else
                        trialClock{j,k}(i)=NaN;
                    end
                    
                    %Previous is counterclockwise to current then to prior
                    %--------------------------------------------------
                    %get current trials when the sequence is: prior < current
                    %<previous.We set -179 and not -180 to be sure they are not
                    %at opposite direction in which case interpretation is not
                    %possible.
                    if  angle.PriortoCur(i)<0 && angle.PriortoCur(i)>-179 && angle.CurtoPrev(i)<0 && angle.CurtoPrev(i)>-179
                        
                        %store the position of those trials
                        trialCClock{j,k}(i)=posCurrentThisCondition(i);
                    else
                        trialCClock{j,k}(i)=NaN;
                    end
                    
                end
                
                %get rid of NaN
                trialClock{j,k}(isnan(trialClock{j,k}))=[];
                trialCClock{j,k}(isnan(trialCClock{j,k}))=[];
            end
        end
    end
end

%(case bimodal prior)
%--------------------
if sum(strcmp(varargin{1},'experiment'))
    if strcmp(varargin{1}(find(strcmp(varargin{1},'experiment'))+1),'bimodalPrior')
        %counterclock mode at each trial
        %when the distance between the first and second mode is <0 the
        %first mode is more counterclockwise than the second and when it is
        %>0, the second mode is more clockwise than the first. Now get the
        %coordinates of the counterclock mode (left)
        distModes=SLvectors2signedAngle(FA.g2.thisT(:,1),FA.g2.thisT(:,2));
        cclockMode(distModes<0)=FA.g2.thisT(distModes<0,1);
        cclockMode(distModes>0)=FA.g2.thisT(distModes>0,2);
        dir.priorcclockModeCoor=polar2cartesian(cclockMode,1);
        
        %Now get the coordinates of the clock mode (right)
        clockMode(distModes<0)=FA.g2.thisT(distModes<0,2);
        clockMode(distModes>0)=FA.g2.thisT(distModes>0,1);
        dir.priorclockModeCoor=polar2cartesian(clockMode,1);
        
        %case factor 1 is the prior modes
        if strcmp(FA.g1.nm,'priormodes')
            FA.g1.thisTSeqAna=FA.g1.distance_mode;
            %case not prior mode
        else
            FA.g1.thisTSeqAna=FA.g1.thisT;
        end
        
        %case factor 2 is the prior modes
        if strcmp(FA.g2.nm,'priormodes')
            FA.g2.thisTSeqAna=FA.g2.distance_mode;
            %case not prior mode
        else
            FA.g2.thisTSeqAna=FA.g2.thisT;
        end
        
        %case factor 3 is the prior modes
        if strcmp(FA.g3.nm,'priormodes')
            FA.g3.thisTSeqAna=FA.g3.distance_mode;
            %case not prior mode
        else
            FA.g3.thisTSeqAna=FA.g3.thisT;
        end
        
        %calculate and draw
        for j=1:FA.g1.lvlsnb
            
            %find this factors' trials
            trialsThisFA1=FA.g1.thisTSeqAna==FA.g1.lvlsnm(j);
            
            %look at factor 2
            for k=1:FA.g2.lvlsnb
                
                %find this factor's trials
                trialsThisFA2=FA.g2.thisTSeqAna==FA.g2.lvlsnm(k);
                
                %get positions of current trials with this condition
                posCurrentThisCondition=find(trialsThisFA1 & trialsThisFA2);
                
                %get positions of the trials that preceded
                posPrevious=posCurrentThisCondition-1;
                
                %remove first trials because there is not a preceding trial
                posPrevious(posCurrentThisCondition==1)=[];
                posCurrentThisCondition(posCurrentThisCondition==1)=[];
                
                %get current trial's motion direction
                dir.dispCoorCurThisCondition=dir.dispCoor(posCurrentThisCondition,:);
                
                %get counterclock mode in the current trial
                %get clock mode in the current trial
                dir.priorcclockModeCoorThisCondition=dir.priorcclockModeCoor(posCurrentThisCondition,:);
                dir.priorclockModeCoorThisCondition=dir.priorclockModeCoor(posCurrentThisCondition,:);
                
                %get the previous motion direction
                dir.dispCoorPrev=dir.dispCoor(posPrevious,:);
                
                
                %get the current directions
                for i=1:size(dir.dispCoorCurThisCondition,1)
                    
                    %When current and previous directions were displayed
                    %counterclockwise to the most counterclockwise mode of
                    %the bimodal prior (i.e.,0< angle(left prior mode,current)<180)
                    angle.ccPriortoCur(i)=SLvectors2signedAngle(dir.priorcclockModeCoorThisCondition(i,:),...
                        dir.dispCoorCurThisCondition(i,:));
                    
                    %clockwise mode
                    angle.cPriortoCur(i)=SLvectors2signedAngle(dir.priorclockModeCoorThisCondition(i,:),...
                        dir.dispCoorCurThisCondition(i,:));
                    
                    %angle between current and previous direction
                    angle.CurtoPrev(i)=SLvectors2signedAngle(dir.dispCoorCurThisCondition(i,:),...
                        dir.dispCoorPrev(i,:));
                    
                    %previous is clockwise to current and to most clockwise prior
                    %mode
                    %------------------------------------------------------------
                    %get current trials when the sequence is: past < current < prior
                    %case current counterclockwise to prior mean
                    %case previous counterclockwise to current for this condition
                    if  0<angle.ccPriortoCur(i) && angle.ccPriortoCur(i)<180 && 0<angle.CurtoPrev(i) && angle.CurtoPrev(i)<180
                        
                        %store the position of those trials
                        trialClock{j,k}(i)=posCurrentThisCondition(i);
                    else
                        trialClock{j,k}(i)=NaN;
                    end
                    
                    %previous is counterclockwise to current and most
                    %counterclocwise prior
                    %-----------------------------------------------
                    %get current trials when the sequence is: past < current < prior
                    %case current counterclockwise to prior mean
                    %case previous counterclockwise to current for this condition
                    if  angle.cPriortoCur(i)<0 && angle.cPriortoCur(i)>-179 && angle.CurtoPrev(i)<0 && angle.CurtoPrev(i)>-179
                        
                        %store the position of those trials
                        trialCClock{j,k}(i)=posCurrentThisCondition(i);
                    else
                        trialCClock{j,k}(i)=NaN;
                    end
                    
                end
                
                %now get rid of NaN
                trialClock{j,k}(isnan(trialClock{j,k}))=[];
                trialCClock{j,k}(isnan(trialCClock{j,k}))=[];
                
            end
        end
    end
end

%Supplementary analyses
%----------------------
%%Plot the effect of Factor on the data (e.g., prior's strength)
function [fig3] = graphPriorEffect(fig1,FA)
%Calculate the effect of the first FA on the data (level 2 - level 1)
%note: it is possible two compare the data for two levels only.

%set figure
%fig3.hdle=figure('Position', [0 0 1000 400]);%pixels
fig3.hdle=figure('color','w'); %pixels
axis square
fig3.nm=[fig1.nm,'_PriorEffect'];

%Check if the first FA has exactly 2 levels.
if FA.g1.lvlsnb ~=2
    fprintf ('%s \n','(analysis) graphPriorEffect analysis has been cancelled. This analysis is valid only if the number of levels in the first FA=2')
    return
end
fig3.datamean=[fig1.datamean.deg(:,:,1) - fig1.datamean.deg(:,:,2)]'; %e.g., prior's width effect


%Draw linear fit to the data
hold all;
for j=1:FA.g2.lvlsnb
    linefit(FA.g3.lvlsnm,fig3.datamean(:,j),1:10:360,FA.g2.color{j});
end

%Draw the prediction for an absence of effect of the first FA (e.g., prior width)
plot(FA.g3.lvlsnm,zeros(FA.g3.lvlsnb,1),...
    '--','color',[0.5 0.5 0.5]);

%Draw the mean of the prior
plot([140 140],[min(min(fig3.datamean)) max(max(fig3.datamean))],'b:',...
    'linewidth', 1);

%Draw the effect of the first FA(e.g., prior's width)
for j=1:FA.g2.lvlsnb
    p13(j)=plot(FA.g3.lvlsnm,fig3.datamean(:,j),...
        'o',...
        'markersize',12  ,... %18 for illustrator (9 for publishing)
        'markerfacecolor', FA.g2.color{j},...
        'MarkerEdgeColor','w',...
        'color', FA.g2.color{j},...
        'displayname',strcat(FA.g2.nm,'--',num2str(FA.g2.lvlsnm(j))));
end

%Set the legend
leg=legend (p13,'location','northwest');
set(leg, 'Box', 'off');

%Set the axes label
ylabel('Difference of estimated directions (?, narrow - wide)','fontsize',11);
xlabel({'Displayed direction','(degrees)'},'fontsize',11);

%set the graph parameters
%set the unit of the x-axis
xunit=1:6:FA.g3.lvlsnb;
set(gca,...
    'xtick',FA.g3.lvlsnm(xunit),'xticklabel',FA.g3.lvlsnm(xunit),...
    'fontsize',11);

xlim([min(FA.g3.lvlsnm)-10 max(FA.g3.lvlsnm)+10]);

%backup figures and parameters
autobackup(fig3.hdle,fig3.nm,'.fig');




%%Support functions
%------------------
%Backup
function autobackup(task,filename,filetype)
%input:
%structure called 'fig' with 2 fields:
%- 'name': e.g., '120910_Steeve_exp02_resul_run8_var1000_mean195_StimStrength024_t250_sess2'
%- 'hdle': the handle of the figure you want to backup

%check if the file name exists
r=dir;
clear i
nametocheck=strcat(filename,filetype);
for i=1:length(r); scanres(i)=strcmp(r(i).name,nametocheck); end

%if the file name exists, increment
i=0;
if ~isempty(find(scanres==1))
    while ~isempty(find(scanres==1))
        i=i+1;
        %increment the name
        filename=strcat(filename,'_0',num2str(i));
        nametocheck=strcat(filename,filetype);
        %check if the name exists already
        for j=1:length(r); scanres(j)=strcmp(r(j).name,nametocheck); end
    end
    errorms=[' "This filename exists already. The new filename is "',filename,'" '];
    disp(errorms)
end
if strcmp(filetype,'.fig')
    saveas(task,filename,'fig'); %name is a string
elseif strcmp(filetype,'.mat')
    save(f0ilename,inputname(1)); %name is a string
end

%Do analysis of covariance (ancova)
function [Stat] = statAncova(x,y,group)
%
%%Input:
%%x: e.g.,displayed dir
%%y: e.g.,estimated dir.
%%group: "strong prior" vs "weak prior".
%
%%OUTPUT: interaction p-value in Anova table specifies if slopes differ.
%
%%Important issue! I can't find how to apply ancova on circular data.
%%Solution: Ancova can be performed with the vector averaged displayed estimated
%%dirs and the displayed dirs.
%
%%x
%databank.data(:,8)
%
%%y
%
%%group
%databank.data(:,6)
%
%check enough group
if length(unique(group))<2
    Stat = [];
    sprintf('(statAncova) Not enough groups for ANCOVA')
    return
else
    [h,atab,ctab,stats] = aoctool(x,y,group,0.05,'','','','off');
    Stat={'h','anovatab','coeftab','stats'};
    Stat(2,:)={h,atab,ctab,stats};
end

%Draw the raw data (e.g., individual dir estimates)
function [fig4] = drawRawEstDir(fig2)

%still under developpment

%input
%angles=vector of cartesian coordinates on the unit circle.
%r=radius of the motion patch: 2.5 degrees.

%Title the panels
figTitles={'StimStrength24','StimStrength12','StimStrength6'};

%set the factors
%factor 1(e.g., Stim strength )
FA1.lvlnum=size(fig2.datacoord, 1);
%factor 2(e.g., displayed dir)
FA2.lvlnum=size(fig2.datacoord, 2);
%factor 3(e.g., prior)
FA3.lvlnum=size(fig2.datacoord, 3);

%Initialize the figure
%loop over factor 2 (e.g., coherence)
for  k=1:FA3.lvlnum
    fig4(k)=figure('color','w');
    %loop over factor 3 (e.g., displayed dirs)
    for i=1:FA2.lvlnum
        r=2.5;
        %initialize figure
        subplot(6,6,i)
        
        %intitialize polar
        p=polar(0,r);
        set(p, 'Visible', 'Off');
        
        %%    Remove unnecessary lines and text.
        %    %Find the lines in the polar plot
        %    h=findall(gcf,'type','line');
        %    %Remove the handle for the polar plot line from the array
        %    h(h ==p)=[];
        %%    Delete other lines
        %    delete(h);
        %    %Find all text objects in the polar plot
        %    t=findall(gcf,'type','text');
        %    %Delete text objects
        %    delete(t);
        
        %set arrows
        arrowsz=2;
        arrowidth=1;
        
        %        FA1.lvlnum=size(figOut{2}.datacoord, 1);
        %        FA2.lvlnum=size(figOut{2}.datacoord, 2);
        %        FA3.lvlnum=size(figOut{2}.datacoord, 3);
        
        %set the colors of factor 3
        c=colormap;
        FA3.lvlcolor={[0.5 0 0],...
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
        
        %set parameters for the graphics
        a=[1.3; 1.6; 1.9; 2.1; 2.4; 2.7];
        
        %loop over the levels of factor 3 (e.g., stim strength)
        for j=1:FA1.lvlnum
            cartCoord.all =fig2.datacoord(j,i,k).all;
            cartCoord.mean=fig2.datacoord(j,i,k).mean;
            
            arrow_lengthPrev=sqrt(cartCoord.mean(1)^2 + cartCoord.mean(2)^2);
            scale2radius=2.5/arrow_lengthPrev;
            
            %draw raw and mean estimated dirs
            if ~isempty(cartCoord.all)
                arrow(0.9*cartCoord.all*a(j),cartCoord.all*a(j), arrowsz, [], [], arrowidth,'EdgeColor', FA3.lvlcolor{j},'facecolor', FA3.lvlcolor{j})
                arrow([0,0],scale2radius*cartCoord.mean, 2, [], [], 1,'EdgeColor', FA3.lvlcolor{j},'facecolor', FA3.lvlcolor{j})
            end
            %draw displayed dirs
            hold on;
            h2=polar(fig2.dataInfog3{j,i,k}*pi/180, 2*r, 'gs');
            set(h2,'markerfacecolor',[0 0.65 0], 'Markersize',10)
            %            leg=legend (p12,'location','northwest');
            %            set(leg, 'Box', 'off');
        end
    end
    title (figTitles{k})
    fprintf('(analyses) Done ! \n')
    
end




