
%SLanalysesEyeMvt.m
%
%
% author: steeve laquitaine
%   date: 140911 updated 150126
%purpose: analyse eye movements
%
%usage:
%
%ex1:
%
%       datapath = slgetdatapath('lab')      
%       [eyedatabank,output] = SLanalysesEyeMvt({'sub01'},...
%                   'dataPath',datapath,...
%                   'AnovaEyeMvt',...
%                   'plotEyeMvtStats',....
%                   'filename','test')
% 
%
%ex2:
%
%       %Quick check one stim file
%       %Move one stim file (.mat and .edf) to ~/desktop 
%       SLanalysesEyeMvt
%
%
%
%list of analyses
%----------------
%   - anova of eye movement per conditions ('AnovaEyeMvt')
%   - plot of eye movement per conditions ('plotEyeMvtStats')
%
%
%varargin
%--------
%   select a session
%       'session',5,...
%
%   anova of eye movement per conditions
%       'AnovaEyeMvt'
%
%   plot of eye movement per conditions
%       'plotEyeMvtStats'
%
%
%Requires mgl distribution (Justin Gardner)
%------------------------------------------
%       getTaskEyeTraces
%
%To do:
%------
%  We should probably correct with bonferroni multiple comparison because we
%  have multiple factors (4 priors, 12 comparisons);
%
%  Quick load stimfile (e.g.,stim01.edf) from desktop.

function [eyedatabank,output] = SLanalysesEyeMvt(subjects,varargin)

%------------------- Quick check one stim file ------------------
%data are on Desktop
if nargin==0
    slPrintfStr('analyses',' Checking data...')
    %get data
    cd ~/Desktop/   
    %default analyses
    varargin = {'AnovaEyeMvt','plotEyeMvtStats','filename','eyedatbank'};
       
    %loop over the files and collect their data
    eyedatabank.pstdForTraj = nan(1,1);
    
    %edf filename
    datalist = dir('*.edf');
    if isempty(datalist)
        slPrintfStr('SLanalysesEyeMvt',' No .edf found on ~/Desktop')
        keyboard
    elseif length(datalist) > 1
        slPrintfStr('SLanalysesEyeMvt',' More than 1 .edf found on ~/Desktop')
        keyboard
    end
    datalisting{1} = datalist.name;
    thisfilename{1} = datalisting{1}(1:end-4);
    
    %eye data starting from motion onset
    %requires stim..mat to be in ~/Desktop
    %i.e., segNum=2
    eye_data{1} = SLgetTaskEyeTraces(thisfilename{1},...
        'myRandomDir_x_myRandomCoh',...
        'taskNum=2','removeBlink=1','EdfRenamed',...
        'dispFig=0','segNum=2');
    
    %period eye data fixation start to before end of motion.
    %i.e., segNum=1 and time <1.3s
    %eyeTend = find(eye_data{i}.eye.time < 1.3);
    eyeTend = find(eye_data{1}.eye.time < 0.3);
    eyexPos{1} = nanmean(eye_data{1}.eye.xPos(:,eyeTend),2);
    eyeyPos{1} = nanmean(eye_data{1}.eye.yPos(:,eyeTend),2);
    eyetime{1} = eye_data{1}.eye.time(eyeTend);
    
    %eye trajectories
    eyedatabank.xTraj{1} = eye_data{1}.eye.xPos(:,eyeTend);
    eyedatabank.yTraj{1} = eye_data{1}.eye.yPos(:,eyeTend);
    
    %conditions
    %----------
    %direction and coh
    if isfield(eye_data{1}.randVars,'myRandomDir')
        direction{1} = SLmakeColumn(eye_data{1}.randVars.myRandomDir);
    elseif isfield(eye_data{1}.randVars,'myRandomloc')
        direction{1} = SLmakeColumn(eye_data{1}.randVars.myRandomloc);
    end
    if isfield(eye_data{1}.randVars,'myRandomCoh')
        coh{1} = SLmakeColumn(eye_data{1}.randVars.myRandomCoh);
    elseif isfield(eye_data{1}.randVars,'myRandomCon')
        coh{1} = SLmakeColumn(eye_data{1}.randVars.myRandomCon);
    end
    
    %prior std in filename
    filename = thisfilename{1};
    load(thisfilename{1});
    datai = getTaskParameters(myscreen,task);
    
    %get the std of the prior from the file nm
    if isempty(isfield(datai{2}.randVars,'myStrength'))
        fprintf('%s /n','(SLMakedatabank) The filename does not contains "myStrength" information...')
        keyboard
    end        
    if floor(datai{2}.randVars.myStrength*10)/10==0.7
        Pstd_thisT = 80;
    elseif floor(datai{2}.randVars.myStrength*10)/10==2.7
        Pstd_thisT = 40;
    elseif floor(datai{2}.randVars.myStrength*10)/10==8.7
         Pstd_thisT = 20;
    elseif floor(datai{2}.randVars.myStrength*10)/10==33.3
         Pstd_thisT = 10;
    end
    numTrials = length(datai{2}.randVars.myStrength);
    pstd = repmat(Pstd_thisT,numTrials,1);    
    eyedatabank.pstdForTraj(1) =  pstd(1);
    session{1} = repmat(1,numel(direction{1}),1);
    
    %databank
    eyedatabank.d = direction{1};
    eyedatabank.coh = coh{1};
    eyedatabank.pstd = pstd;
    eyedatabank.eyexPos = eyexPos{1};
    eyedatabank.eyeyPos = eyeyPos{1};
    eyedatabank.eyeTime = eyetime;
    eyedatabank.session = session{1};
        
else
    
    %------------------- subject-based analysis---------------------
    
    %call for help
    if ieNotDefined('subjects')
        help analyses
        return
    end
    
    %backup this .m file in the worspace
    %-----------------------------------
    mfilename = SLgetActivemFile;
    myMfile = SLbackupMfileInWSpace(mfilename);
    output.mfilename  = mfilename;
    output.myMfile = myMfile;
    
    %get data path
    if sum(strcmp(varargin,'dataPath'))
        dataPath = varargin{find(strcmp(varargin,'dataPath'))+1};
        cd(dataPath)
    else
        fprintf('%s \n','(SLfitBayesianModel) You need to set dataPath')
        return
    end
    
    %get the name under which we save the file
    if sum(strcmp(varargin,'filename'))
        filename = varargin{find(strcmp(varargin,'filename'))+1};
    else
        fprintf('%s \n',['(SLfitBayesianModel) You need to set the name',...
            ' of the .mat file that will be saved'])
    end
    
    %case we want one session only
    if sum(strcmp(varargin,'session'))
        selSessions = varargin{find(strcmp(varargin,'session'))+1};
    end
    
    %set the directory to analyse
    %case we are already in data directory go to parent
    a=pwd;
    if strcmp(a(end-4:end),'/data')
        pathFold=pwd;
    else
        %otherwise go to /data directory
        pathFold=strcat(pwd,'/data');
    end
    d=dir(pathFold);
    
    %collect subjects
    datatmp.subjects=subjects;
    
    %find the subfolders to analyse
    for ii = 1 : numel(datatmp.subjects)
        isubType(ii) = find(strcmp({d.name},datatmp.subjects(ii)));
    end
    
    %get the names of the subfolders
    Folds.name = {d(isubType).name}';
    
    %check if correct
    disp(Folds.name)
    %uiwait(msgbox('Check out the command window for the folder names.'));
    
    %count the number of subfolders
    Folds.nb = numel(Folds.name);
    
    %initialize databank
    eyedatabank.d = [];
    eyedatabank.coh = [];
    eyedatabank.pstd = [];
    eyedatabank.eyexPos = [];
    eyedatabank.eyeyPos = [];
    eyedatabank.eyeTime = [];
    eyedatabank.session  = [];
    
    %loop over the subjects folders
    for j = 1 : Folds.nb
        
        %directory of the subfolder to open
        Folds.path(j)=strcat(pathFold,'/',Folds.name(j));
        
        %switch to the subfolder to open
        cd(Folds.path{j})
        
        %delete svn files (if exist)
        !rm ._*;
        
        %specify the files to load
        datadir=dir('*data*.mat');%directory
        
        %check if there are data in the directory
        if isempty(datadir)
            disp(strcat('No data were found in directory: ',Folds.name(j)))
            return
        end
        
        %collect the name of each file
        datatmp.nm={datadir.name};
        
        %count the number of file
        datatmp.nb=numel(datatmp.nm);
        
        %check if data have been specified in the function argin
        if ~isempty(datatmp.nm)
            datalisting = datatmp.nm;
        else
            %remove possible svn files
            datalisting = dir('*data*.mat');
            datalisting = {datalisting.name};
        end
        
        %case we want a specific session
        %-------------------------------
        if sum(strcmp(varargin,'session'))
            
            
            %find the files for this session
            posFile =[];
            for ijk = 1 : numel(datalisting)
                
                %files that match our session
                sessiontm = str2double(datalisting{ijk}(strfind(datalisting{ijk},...
                    'sess')+5));
                if sessiontm == selSessions
                    posFile = [posFile ; ijk];
                end
            end
            
            %Status
            fprintf('%s \n','(SLanalysesEyeMvt) Selected sessions are:')
            fprintf('%s \n','---------------------------------------------')
            for iSel = 1: numel(posFile)
                fprintf('%s \n',datalisting{posFile(iSel)})
            end
            fprintf('%s \n','---------------------------------------------')
            
            %data from selected session files
            eyedatabank.pstdForTraj = nan(1,numel(posFile));
            for iSel = 1 : numel(posFile)
                
                %session files
                ThisFile = posFile(iSel);
                
                %edf data
                load(datalisting{ThisFile});
                if isfield(myscreen.eyetracker,'datafilename')
                    
                    %edf filename
                    thisfilename{iSel} = datalisting{ThisFile}(1:end-4);
                    
                    %eye data starting from motion onset
                    %i.e., segNum=2
                    eye_data{iSel} = SLgetTaskEyeTraces(thisfilename{iSel},...
                        'myRandomDir_x_myRandomCoh',...
                        'taskNum=2','removeBlink=1','EdfRenamed',...
                        'dispFig=0','segNum=2');
                    
                    if isfield(eye_data{iSel},'eye')
                        
                        %period is during motion stimulus presentation (300 ms)
                        %we take the mean position of the eye during motion period
                        %for each trial (We later plot the mean of this mean
                        %position over trials and the variability over trials).
                        eyeTend = find(eye_data{iSel}.eye.time < 0.3);
                        eyexPos{iSel} = nanmean(eye_data{iSel}.eye.xPos(:,eyeTend),2);
                        eyeyPos{iSel} = nanmean(eye_data{iSel}.eye.yPos(:,eyeTend),2);
                        eyetime{iSel} = eye_data{iSel}.eye.time(eyeTend);
                        
                        %eye trajectories
                        eyedatabank.xTraj{iSel} = eye_data{iSel}.eye.xPos(:,eyeTend);
                        eyedatabank.yTraj{iSel} = eye_data{iSel}.eye.yPos(:,eyeTend);
                        
                        %conditions
                        %----------
                        %direction and coh
                        direction{iSel} = SLmakeColumn(eye_data{iSel}.randVars.myRandomDir);
                        coh{iSel} = SLmakeColumn(eye_data{iSel}.randVars.myRandomCoh);
                        
                        %prior std in filename
                        if isempty(strfind(thisfilename{iSel},'Pstd'))
                            fprintf('%s /n','(SLanalysesEyeMvt) the filename does not contains "Pstd" information...')
                            return
                        end
                        pstdtmp = str2double(thisfilename{iSel}(strfind(thisfilename{iSel},...
                            'Pstd')+numel('Pstd'):strfind(thisfilename{iSel},...
                            'Pstd')+numel('Pstd')+2));
                        pstd{iSel} = repmat(pstdtmp,numel(direction{iSel}),1);
                        eyedatabank.pstdForTraj(iSel) = pstdtmp;
                        
                        %get the session
                        if isempty(strfind(datalisting{ThisFile},'sess'))
                            fprintf('%s /n','(SLMakedatabank) the filename does not contains "session" information...')
                            return
                        end
                        sessiontm = str2double(datalisting{ThisFile}(strfind(datalisting{ThisFile},...
                            'sess')+5));
                        session{iSel} = repmat(sessiontm,numel(direction{iSel}),1);
                    else
                        %empty
                        origStimFileName{ijk}  = [];
                        eyexPos{ijk} = [];
                        eyeyPos{ijk} = [];
                        eyetime{ijk} = [];
                        direction{ijk} = [];
                        coh{ijk} = [];
                        pstd{ijk} = [];
                        session{ijk} = [];
                    end
                else
                    %empty
                    origStimFileName{ijk}  = [];
                    eyexPos{ijk} = [];
                    eyeyPos{ijk} = [];
                    eyetime{ijk} = [];
                    direction{ijk} = [];
                    coh{ijk} = [];
                    pstd{ijk} = [];
                    session{ijk} = [];
                end
                
                clear myscreen
                
                %databank
                eyedatabank.d = [eyedatabank.d ; direction{iSel}];
                eyedatabank.coh = [eyedatabank.coh ; coh{iSel}];
                eyedatabank.pstd = [eyedatabank.pstd ; pstd{iSel}];
                eyedatabank.eyexPos = [eyedatabank.eyexPos ; eyexPos{iSel}];
                eyedatabank.eyeyPos = [eyedatabank.eyeyPos ; eyeyPos{iSel}];
                eyedatabank.eyeTime = eyetime;
                eyedatabank.session = [eyedatabank.session ; session{iSel}];
            end
            
        else
            
            %case we want all sessions (default)
            %----------------------------------
            %loop over the files and collect their data
            eyedatabank.pstdForTraj = nan(1,numel(datalisting));
            
            for ijk = 1 : numel(datalisting)
                
                %edf data
                load(datalisting{ijk});
                if isfield(myscreen.eyetracker,'datafilename')
                    
                    %edf filename
                    thisfilename{ijk} = datalisting{ijk}(1:end-4);
                    
                    %eye data starting from motion onset
                    %i.e., segNum=2
                    eye_data{ijk} = SLgetTaskEyeTraces(thisfilename{ijk},...
                        'myRandomDir_x_myRandomCoh',...
                        'taskNum=2','removeBlink=1','EdfRenamed',...
                        'dispFig=0','segNum=2');
                    
                    if isfield(eye_data{ijk},'eye')
                        
                        %period eye data fixation start to before end of motion.
                        %i.e., segNum=1 and time <1.3s
                        %eyeTend = find(eye_data{i}.eye.time < 1.3);
                        eyeTend = find(eye_data{ijk}.eye.time < 0.3);
                        eyexPos{ijk} = nanmean(eye_data{ijk}.eye.xPos(:,eyeTend),2);
                        eyeyPos{ijk} = nanmean(eye_data{ijk}.eye.yPos(:,eyeTend),2);
                        eyetime{ijk} = eye_data{ijk}.eye.time(eyeTend);
                        
                        %eye trajectories
                        eyedatabank.xTraj{ijk} = eye_data{ijk}.eye.xPos(:,eyeTend);
                        eyedatabank.yTraj{ijk} = eye_data{ijk}.eye.yPos(:,eyeTend);
                        
                        %conditions
                        %----------
                        %direction and coh
                        direction{ijk} = SLmakeColumn(eye_data{ijk}.randVars.myRandomDir);
                        coh{ijk} = SLmakeColumn(eye_data{ijk}.randVars.myRandomCoh);
                        
                        %prior std in filename
                        if isempty(strfind(thisfilename{ijk},'Pstd'))
                            fprintf('%s /n','(SLanalysesEyeMvt) the filename does not contains "Pstd" information...')
                            return
                        end
                        pstdtmp = str2double(thisfilename{ijk}(strfind(thisfilename{ijk},...
                            'Pstd')+numel('Pstd'):strfind(thisfilename{ijk},...
                            'Pstd')+numel('Pstd')+2));
                        pstd{ijk} = repmat(pstdtmp,numel(direction{ijk}),1);
                        eyedatabank.pstdForTraj(ijk) = pstdtmp;
                        
                        %get the session
                        if isempty(strfind(datalisting{ijk},'sess'))
                            fprintf('%s /n','(SLMakedatabank) the filename does not contains "session" information...')
                            return
                        end
                        sessiontm = str2double(datalisting{ijk}(strfind(datalisting{ijk},...
                            'sess')+5));
                        session{ijk} = repmat(sessiontm,numel(direction{ijk}),1);
                    else
                        %empty
                        origStimFileName{ijk}  = [];
                        eyexPos{ijk} = [];
                        eyeyPos{ijk} = [];
                        eyetime{ijk} = [];
                        direction{ijk} = [];
                        coh{ijk} = [];
                        pstd{ijk} = [];
                        session{ijk} = [];
                    end
                else
                    %empty
                    origStimFileName{ijk}  = [];
                    eyexPos{ijk} = [];
                    eyeyPos{ijk} = [];
                    eyetime{ijk} = [];
                    direction{ijk} = [];
                    coh{ijk} = [];
                    pstd{ijk} = [];
                    session{ijk} = [];
                end
                clear myscreen
                
                %databank
                eyedatabank.d = [eyedatabank.d ; direction{ijk}];
                eyedatabank.coh = [eyedatabank.coh ; coh{ijk}];
                eyedatabank.pstd = [eyedatabank.pstd ; pstd{ijk}];
                eyedatabank.eyexPos = [eyedatabank.eyexPos ; eyexPos{ijk}];
                eyedatabank.eyeyPos = [eyedatabank.eyeyPos ; eyeyPos{ijk}];
                eyedatabank.eyeTime = eyetime;
                eyedatabank.session = [eyedatabank.session ; session{ijk}];                
            end
        end
    end
end


%save file in 'eyeAnal' folder
mkdir('eyeAnal')
cd eyeAnal
    
%case Anova of eye movement
output=[];
if sum(strcmp(varargin,'AnovaEyeMvt'))
    output = AnovaEyeMvt(eyedatabank,output,'prior');
end

%case Anova of eye movement
if sum(strcmp(varargin,'plotEyeMvtStats'))
    output = plotEyeMvtStats(eyedatabank,output,'prior');
    %save image
    %print('-depsc','-r300',filename)
end

%save file
save(filename)



%plot eye movements
function output = plotEyeMvtStats(eyedatabank,output,VarName)

%prior condition
if strcmp(VarName,'prior')
    
    figure('color','w')
    set(gcf,'position',[680 575 507 523])
    colors = [0.5 0 0;...
        1 0.0625 0;...
        1 0.4 0;...
        0 0.7 0];
    
    %elliptic space
    kk = -pi:0.01:pi;
    
    %priors
    pstds = unique(eyedatabank.pstd);
    pstds = sort(pstds,'descend');
    
    for i = 1 : numel(pstds)
        thisPstd = find(eyedatabank.pstd==pstds(i));
        
        %means position across trials
        output.meanEyexPos(i)= nanmean(eyedatabank.eyexPos(thisPstd));
        output.meanEyeyPos(i)= nanmean(eyedatabank.eyeyPos(thisPstd));
        
        %std
        output.stdEyexPos(i)= nanstd(eyedatabank.eyexPos(thisPstd));
        output.stdEyeyPos(i)= nanstd(eyedatabank.eyeyPos(thisPstd));
        
        %ellipses
        x(:,i) = output.meanEyexPos(i) + output.stdEyexPos(i)*cos(kk);
        y(:,i) = output.meanEyeyPos(i) + output.stdEyeyPos(i)*sin(kk);
        
        hold all
        
        %mean positions and variability
        %------------------------------
        %area
        plot(x(:,i),y(:,i),'color',colors(i,:),'linesmoothing','on')
        
        %mean
        plot(output.meanEyexPos(i),output.meanEyeyPos(i),'.',...
            'color',colors(i,:),'markersize',25);
    end
    %legend
    legend('80 deg pstd','mean','40 deg pstd','mean','20 deg pstd','mean',...
        '10 deg pstd','mean')
    legend('boxoff')
    
    %fixation axis
    plot(linspace(0,0,100),linspace(-1.5,1.5,100),'k:')
    plot(linspace(-1.5,1.5,100),linspace(0,0,100),'k:')
    
    %stimulus borders
    plot(2.5*cos(kk),2.5*sin(kk),'k:','linesmoothing','on')
    
    %graphics
    %radius stimulus is 2.5 deg
    %xlim([-2.5 2.5])
    %ylim([-2.5 2.5])
    xlabel('horizontal movements (deg)')
    ylabel('vertical movements (deg)')
    
    %Anova of conditions
    %case stat have not been calculated
    if isempty(output)
        output = AnovaEyeMvt(eyedatabank,output,VarName);
    end
    
    %title and info
    sessions = num2str(unique(eyedatabank.session));
    dfxGroup = num2str(output.tabx{2,3});
    dfxWithin = num2str(output.tabx{3,3});
    Fvalx = num2str(output.tabx{2,5});
    dfyGroup = num2str(output.taby{2,3});
    dfyWithin = num2str(output.taby{3,3});
    Fvaly = num2str(output.taby{2,5});
    
    title({['Session ' sessions'],...
        'Anova(''priors'') for x and y coord',...
        ['x: [F(',dfxGroup,',',dfxWithin,')=',Fvalx,', ',...
        'p=',sprintf('%.1g ',output.Pvalx) ']'],...
        ['y: [F(',dfyGroup,',',dfyWithin,')=',Fvaly,', ',...
        'p=',sprintf('%.1g ',output.Pvaly) ']']})
    
    text(0,-2.5,'Stimulus border')
    axis square
end

%Anova of conditions
function output = AnovaEyeMvt(eyedatabank,output,VarName)

%prior condition
if strcmp(VarName,'prior')
    
    %p-values
    [output.Pvalx,output.tabx,output.statx] = anova1(eyedatabank.eyexPos,eyedatabank.pstd);
    [output.Pvaly,output.taby,output.statY] = anova1(eyedatabank.eyeyPos,eyedatabank.pstd);
end

%every combination of condition
if strcmp(VarName,'_every_')
    
    %conditions and their trial position
    [cond,b,c] = SLuniqpair([eyedatabank.d eyedatabank.coh eyedatabank.pstd]);
    [~,pos]= sortrows(c);
    
    %sort data per condition
    conditionSorted = c(pos);
    eyexPosSorted = eyedatabank.eyexPos(pos);
    eyeyPosSorted = eyedatabank.eyeyPos(pos);
    
    %p-values
    output.Pvalx = anova1(eyexPosSorted,conditionSorted);
    output.Pvaly = anova1(eyeyPosSorted,conditionSorted);
end