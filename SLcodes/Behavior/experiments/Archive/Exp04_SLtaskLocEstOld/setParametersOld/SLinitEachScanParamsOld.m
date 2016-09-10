
%SLinitEachScanParams.m
%
%  author: steeve laquitaine
%    date: 2015/04/23
% purpose: initialize the tmp parameters of several scans such that
%          the prior distribution over locations is spread over all the scans.
%
%          The ouputs "task.mat" are used as input to SLtaskLocEst.m
%
%       usage:
%
% 6 directions 
%              task = SLinitEachScanParams(45:60:355,225,0.74848,[0.06 1],[.5 0 0],...
%                                          'ScanSet1',{'numMain_con1',9,'numMain_con2',6,'numCatch_con1',9,'numCatch_con2',6,'nScans',4},...
%                                          'ScanSet2',{'numMain_con1',7,'numMain_con2',5,'numCatch_con1',6,'numCatch_con2',5,'nScans',1});
% 
%              task = SLinitEachScanParams(45:60:355,225,2.7714,[0.06 1],[1 0.6 0],...
%                                          'ScanSet1',{'numMain_con1',9,'numMain_con2',6,'numCatch_con1',9,'numCatch_con2',6,'nScans',4},...
%                                          'ScanSet2',{'numMain_con1',7,'numMain_con2',5,'numCatch_con1',6,'numCatch_con2',5,'nScans',1});
% 
% 9 directions
%              task = SLinitEachScanParams(25:40:355,225,0.74848,[0.06 1],[.5 0 0],...
%                                          'ScanSet1',{'numMain_con1',9,'numMain_con2',6,'numCatch_con1',9,'numCatch_con2',6,'nScans',4},...
%                                          'ScanSet2',{'numMain_con1',9,'numMain_con2',4,'numCatch_con1',7,'numCatch_con2',4,'nScans',1});
%
%              task = SLinitEachScanParams(25:40:355,225,2.7714,[0.06 1],[1 0.6 0],...
%                                          'ScanSet1',{'numMain_con1',9,'numMain_con2',6,'numCatch_con1',9,'numCatch_con2',6,'nScans',4},...
%                                          'ScanSet2',{'numMain_con1',7,'numMain_con2',5,'numCatch_con1',7,'numCatch_con2',5,'nScans',1});
%
% 18 directions
%              task = SLinitEachScanParams(5:20:355,225,0.74848,[0.06 1],[.5 0 0],...
%                                          'ScanSet1',{'numMain_con1',9,'numMain_con2',6,'numCatch_con1',8,'numCatch_con2',6,'nScans',4},...
%                                          'ScanSet2',{'numMain_con1',6,'numMain_con2',6,'numCatch_con1',9,'numCatch_con2',6,'nScans',1});
%
%              task = SLinitEachScanParams(5:20:355,225,2.7714,[0.06 1],[1 0.6 0],...
%                                          'ScanSet1',{'numMain_con1',9,'numMain_con2',6,'numCatch_con1',8,'numCatch_con2',6,'nScans',4},...
%                                          'ScanSet2',{'numMain_con1',6,'numMain_con2',4,'numCatch_con1',9,'numCatch_con2',4,'nScans',1});

function task = SLinitEachScanParams(samples,Pmean,Pk,cont,colors,varargin)

%move to TaskParameters directory
myPath = mfilename('fullpath');
cd(fileparts(myPath))
cd ../TaskParameters/

%check if task parameters already
%exist in the directory
d = dir('taskScan*');
myPrior = num2str(Pk);
myPrior(myPrior=='.') = '_';

exist = 0;
for i = 1: length(d)
    
    %if already exist, replace or not
    match = SLtextIsContained(d(i).name,['Prior',myPrior]);
    if match == 1
        exist = 1;
    end
end

if exist == 1
    x = input('Scan parameters already exist. Do you want to replace them (y/n)?','s');
    
    %replace existing parameter files
    if x=='y'
       
        delete('taskScan*')
        
        %or abort
    elseif x=='n'
        task = [];
        fprintf('Aborted \n')
        return
    end
end


%generate catch and main trial types
fprintf('(SLinitRunUniPriorLocTask) I am generating parameters for "Main" and "Catch" trials \n')

tmpCatch = SLinitRunUniPriorLocTask('Catch',Pmean,Pk,cont,[43 29],samples,colors); 
SLgraphBigTitle('Prior for Catch Trials, over all scans')

tmpMain = SLinitRunUniPriorLocTask('Main',Pmean,Pk,cont,[44 29],samples,colors); 
SLgraphBigTitle('Prior for Main Trials, over all scans')

delete('Catch*')
delete('Main*')

%Concatenate all conditions within prior
numMain = numel(tmpMain.parameter.loc.series);
numCatch = numel(tmpCatch.parameter.loc.series);

tmp.parameter.loc.trialType = [ones(1,numMain) zeros(1,numCatch)];
tmp.parameter.loc.uniqtrialType = fliplr(unique(tmp.parameter.loc.trialType));
tmp.parameter.loc.ntrialType= numel(tmp.parameter.loc.uniqtrialType);
tmp.parameter.loc.series    = [tmpMain.parameter.loc.series tmpCatch.parameter.loc.series];
tmp.parameter.loc.con       = [tmpMain.parameter.loc.con    tmpCatch.parameter.loc.con];

%Prior
%-----
%check trial types are drawn from the same prior
if tmpMain.parameter.loc.std == tmpCatch.parameter.loc.std
    tmp.parameter.loc.strength  = tmpMain.parameter.loc.std;
else
    fprintf('(SLinitRunUniPriorLocTask) ".strength" differ between Main and Catch trials. The two trial types must be drawn from the same prior... \n')
end

tmp.parameter.loc.inputTrialnumperCon = tmpMain.parameter.loc.inputTrialnumperCon + tmpCatch.parameter.loc.inputTrialnumperCon;
tmp.parameter.loc.TrueTrialnumperCon  = tmpMain.parameter.loc.TrueTrialnumperCon  + tmpCatch.parameter.loc.TrueTrialnumperCon;
tmp.parameter.loc.trueStdperCon = [tmpMain.parameter.loc.trueStdperCon; tmpCatch.parameter.loc.trueStdperCon];

for i = 1: 2
    tmp.parameter.loc.countperCon(i,:) = tmpMain.parameter.loc.countperCon(i,:) + tmpCatch.parameter.loc.countperCon(i,:);
end

tmp.parameter.loc.pperCon         = [{tmpMain.parameter.loc.pperCon} {tmpCatch.parameter.loc.pperCon}];
tmp.parameter.loc.contDistInCount = [{tmpMain.parameter.loc.contDistInCount} {tmpCatch.parameter.loc.contDistInCount}];

%check that main and catch 
%have same mean
if tmpMain.parameter.loc.mean == tmpCatch.parameter.loc.mean
    tmp.parameter.loc.mean = tmpMain.parameter.loc.mean;
else
    fprintf('(SLinitRunUniPriorLocTask) ".mean" differ between Main and Catch trials. The two trial types must be drawn from the same prior... \n')
end

%check that main and catch 
%have same modes
if isfield(tmpMain.parameter.loc,'modes') && isfield(tmpCatch.parameter.loc,'modes')
    if tmpMain.parameter.loc.mean == tmpCatch.parameter.loc.mean
        
        tmp.task.parameter.loc.modes = tmpMain.task.parameter.loc.modes;
        
    end 
else
    fprintf('(SLinitRunUniPriorLocTask) "I couldnt find the prior mode. It will not be saved. \n')
end

%check that main and catch 
%have same samples
if SLisSame(tmpMain.parameter.loc.sample.degree,tmpCatch.parameter.loc.sample.degree)
    tmp.parameter.loc.sample     = tmpMain.parameter.loc.sample;
    tmp.parameter.loc.sampsiz     = tmpMain.parameter.loc.sampsiz;
else
    fprintf('(SLinitRunUniPriorLocTask) ".sample" differ between Main and Catch trials. The two trial types must be drawn from the same prior... \n')
end

tmp.parameter.loc.ntrial = tmpMain.parameter.loc.trialnum + tmpCatch.parameter.loc.trialnum;
tmp.parameter.loc.count  = [{tmpMain.parameter.loc.count} {tmpCatch.parameter.loc.count}];


%Contrast
%--------
%check main and catch have same contrasts
if SLisSame(tmpMain.parameter.loc.conSample,tmpCatch.parameter.loc.conSample) 
    tmp.parameter.loc.uniqCon = tmpMain.parameter.loc.conSample;
else
    fprintf('(SLinitRunUniPriorLocTask) Contrasts differ between Main and Catch. \n')
end
tmp.parameter.loc.nCon = numel(tmp.parameter.loc.uniqCon);


%Organize trials into two scan sets
%with different trial number
%-----------------------------------
if sum(strcmp(varargin,'ScanSet1'))
    
    posParams = find(strcmp(varargin,'ScanSet1'))+1;
    paramSet1 = varargin{posParams};
    
    %Check that all conditions 
    %(trial number) are there
    isexist(1) = sum(strcmp(paramSet1,'numMain_con1'));
    isexist(2)= sum(strcmp(paramSet1,'numMain_con2'));
    isexist(3)= sum(strcmp(paramSet1,'numCatch_con1'));
    isexist(4)= sum(strcmp(paramSet1,'numCatch_con2'));
    
    if sum(isexist)==4
    else
        missingPara = find(isexist~=1);
        fprintf(['The ', num2str(missingPara),' is missing or you may need to check its spelling. \n'])
    end
    
    %get parameter values
    scanSet1.numCond(1) = paramSet1{find(strcmp(paramSet1,'numMain_con1'))+1};
    scanSet1.numCond(2) = paramSet1{find(strcmp(paramSet1,'numMain_con2'))+1};
    scanSet1.numCond(3) = paramSet1{find(strcmp(paramSet1,'numCatch_con1'))+1};
    scanSet1.numCond(4) = paramSet1{find(strcmp(paramSet1,'numCatch_con2'))+1};
end

if sum(strcmp(varargin,'ScanSet1'))
    
    posParams = find(strcmp(varargin,'ScanSet2'))+1;
    paramSet2 = varargin{posParams};
    
    %Check that all conditions 
    %(trial number) are there
    isexist(1) = sum(strcmp(paramSet2,'numMain_con1'));
    isexist(2)= sum(strcmp(paramSet2,'numMain_con2'));
    isexist(3)= sum(strcmp(paramSet2,'numCatch_con1'));
    isexist(4)= sum(strcmp(paramSet2,'numCatch_con2'));
    
    if sum(isexist)==4
    else
        missingPara = find(isexist~=1);
        fprintf(['The ', num2str(missingPara),' is missing or you may need to check its spelling. \n'])
    end
    
    %get trial number per conditon
    scanSet2.numCon(1) = paramSet2{find(strcmp(paramSet2,'numMain_con1'))+1};
    scanSet2.numCon(2) = paramSet2{find(strcmp(paramSet2,'numMain_con2'))+1};
    scanSet2.numCon(3) = paramSet2{find(strcmp(paramSet2,'numCatch_con1'))+1};
    scanSet2.numCon(4) = paramSet2{find(strcmp(paramSet2,'numCatch_con2'))+1};
    
    fprintf('(SLinitRunUniPriorLocTask) You asked for two sets of scans. Parameters differ with between those two sets. \n')
else
    fprintf('(SLinitRunUniPriorLocTask) You asked for one set of scans with the same parameters for all \n')
end

%get number of scans in set 1
if sum(strcmp(varargin,'ScanSet1'))
    posParams = find(strcmp(varargin,'ScanSet1'))+1;
    paramSet1 = varargin{posParams};
    scanSet1.nScans = paramSet1{find(strcmp(paramSet1,'nScans'))+1};
end

%get number of scans in set 2
if sum(strcmp(varargin,'ScanSet2'))
    posParams = find(strcmp(varargin,'ScanSet2'))+1;
    paramSet2 = varargin{posParams};
    scanSet2.nScans = paramSet2{find(strcmp(paramSet2,'nScans'))+1};
end

%check that number of trials 
%in each condition makes sense
ij = 0;
for j = 1 : tmp.parameter.loc.ntrialType
    for i = 1: tmp.parameter.loc.nCon
        ij = ij+1;
        
        %expected
        nTrialsCond(ij) = sum(tmp.parameter.loc.trialType==tmp.parameter.loc.uniqtrialType(j)& ...
            tmp.parameter.loc.con==tmp.parameter.loc.uniqCon(i));
       
    end
end

%input
nInputTrials = scanSet1.nScans * scanSet1.numCond + scanSet2.nScans*scanSet2.numCon;

%Adjust Scan 2 number of trials
if SLisSame(nTrialsCond,nInputTrials)~=1        
        
    scanSet2.numCon = nTrialsCond - scanSet1.nScans * scanSet1.numCond;
   
    fprintf('\n !! WARNING !! The number of trials/condition must be adjusted. Please set scan 2 as follows: !! WARNING !!  \n')
    
    fprintf('------------------------------------------------------------------------------ \n')
    fprintf(['SLinitEachScanParams("ScanSet1",{....},"ScanSet2",{"numMain_con1" ,',num2str(scanSet2.numCon(1)),...
        ' "numMain_con2" ',num2str(scanSet2.numCon(2)),...
        ' "numCatch_con1" ',num2str(scanSet2.numCon(3)),...
        ' "numCatch_con2" ',num2str(scanSet2.numCon(4)),...
        '},"nScans",1) \n']);
    fprintf('------------------------------------------------------------------------------ \n')
    
    keyboard
    
end

%position trials in scan set 1 time series 
tPosBeg = [1 cumsum(scanSet1.numCond(1:end-1))  + 1];
tPosEnd = cumsum(scanSet1.numCond(1:end));

%position trials in scan set 1 time series
tPosBeg2 = [1 cumsum(scanSet2.numCon(1:end-1))  + 1];
tPosEnd2 = cumsum(scanSet2.numCon(1:end));

%count conditions
ij=0;

fprintf('---------------------------------------------    \n')
fprintf('I am allocating conditions to each scan:      \n')
fprintf('---------------------------------------------    \n')

for j = 1 : tmp.parameter.loc.ntrialType
    for i = 1: tmp.parameter.loc.nCon
        
        %count conditions
        ij = ij+1;
        
        %status
        fprintf([' "trialType : ', num2str(tmp.parameter.loc.uniqtrialType(j)), ...
                ' ; Contrast : ', num2str(tmp.parameter.loc.uniqCon(i)),'" \n'])
        
        %get trials this condition
        posCon = find(tmp.parameter.loc.trialType==tmp.parameter.loc.uniqtrialType(j)& ...
            tmp.parameter.loc.con==tmp.parameter.loc.uniqCon(i));
        
        %sample the trials
        posCon = posCon(randperm(numel(posCon)));

        %set trials' position in time series
        tPos = tPosBeg(ij) : tPosEnd(ij); 
        
        %put trials in time series for each scan
        for k = 1 : scanSet1.nScans
            
            %sample trials without replacement
            scanTrials = posCon(1 : scanSet1.numCond(ij));
            
            task{k}.parameter.loc.series(tPos)   = tmp.parameter.loc.series(scanTrials);
            task{k}.parameter.loc.trialType(tPos)= tmp.parameter.loc.trialType(scanTrials);
            task{k}.parameter.loc.con(tPos)      = tmp.parameter.loc.con(scanTrials);
            
            %move to next scan
            posCon(1:scanSet1.numCond(ij)) = [];
        end
        
        
        %Place remaining trials for this 
        %condition in all scans set 2
        
        %status
        fprintf([' "trialType : ', num2str(tmp.parameter.loc.uniqtrialType(j)), ...
                ' ; Contrast : ', num2str(tmp.parameter.loc.uniqCon(i)),'" \n'])
        
        %trials position in time series
        tPos = tPosBeg2(ij) : tPosEnd2(ij);
        
        for l = 1 : scanSet2.nScans
            
            %sample trials without replacement
            if numel(posCon) ~= scanSet2.numCon(ij)
                fprintf(['!! WARNING !! You must adjust the number of trials in ',...
                    'trialType : ', num2str(tmp.parameter.loc.uniqtrialType(j)),...
                    '; Contrast : ',num2str(tmp.parameter.loc.uniqCon(i)),'\n'])
            end
            scanTrials = posCon(1 : scanSet2.numCon(ij));
            
            task{k+l}.parameter.loc.series(tPos)   = tmp.parameter.loc.series(scanTrials);
            task{k+l}.parameter.loc.trialType(tPos)= tmp.parameter.loc.trialType(scanTrials);
            task{k+l}.parameter.loc.con(tPos)      = tmp.parameter.loc.con(scanTrials); 
            
            %move to next scan
            posCon(1:scanSet2.numCon(ij)) = [];
        end 
        
        %status
        if isempty(posCon)
        else
            fprintf(['!!! WARNING !!! Something is wrong. You did not use all the trials...',...
                ' You can still add ',num2str(numel(posCon)),' trials ',...
                'in condition Main:',num2str(tmp.parameter.loc.uniqtrialType(j)),...
                ', cont: ',num2str(tmp.parameter.loc.uniqCon(i)),'\n'])
            keyboard
        end
    end
end
fprintf('done \n')



%fill in other informations
%scan set 1
for k = 1 : scanSet1.nScans
    
    task{k}.parameter.loc.uniqtrialType = unique(task{k}.parameter.loc.trialType);
    task{k}.parameter.loc.strength = tmp.parameter.loc.strength;
    task{k}.parameter.loc.mean = tmp.parameter.loc.mean;
    
    if isfield(tmp.parameter.loc,'modes')
        task{k}.task.parameter.loc.modes = tmp.parameter.loc.modes;
    end
    
    task{k}.parameter.loc.sample = unique(task{k}.parameter.loc.series);
    task{k}.parameter.loc.sampsiz = numel(task{k}.parameter.loc.sample);
    task{k}.parameter.loc.ntrial = numel(task{k}.parameter.loc.series);

end

%scan set 2
for l = 1 : scanSet2.nScans
    
    task{k+l}.parameter.loc.uniqtrialType = unique(task{k+l}.parameter.loc.trialType);
    task{k+l}.parameter.loc.strength = tmp.parameter.loc.strength;
    task{k+l}.parameter.loc.mean = tmp.parameter.loc.mean;
    
    if isfield(tmp.parameter.loc,'modes')
        task{k+l}.parameter.loc.modes = tmp.parameter.loc.modes;
    end
    
    task{k+l}.parameter.loc.sample = unique(task{k+l}.parameter.loc.series);
    task{k+l}.parameter.loc.sampsiz = numel(task{k+l}.parameter.loc.sample);
    task{k+l}.parameter.loc.ntrial = numel(task{k+l}.parameter.loc.series);
end


%graph prior in Main trials for each scan, first contrast
figure('color','w','position',[1236 841 560 259]);
uniqCon = unique(task{1}.parameter.loc.con);
SLgraphBigTitle(['Prior for Main trials, each scan, contrast ',num2str(uniqCon(1))])

nScanTot = scanSet1.nScans + scanSet2.nScans;

for i = 1 : nScanTot
    
    subplot(nScanTot,1,i)
    
    %Main trials and first contrast condition
    MainTrials = task{i}.parameter.loc.trialType==1;
    thisCon = task{i}.parameter.loc.con==uniqCon(1) & MainTrials==1;
    
    hist(task{i}.parameter.loc.series(thisCon),tmpMain.parameter.loc.sample.degree)
    xlim([5 355])
    title(['Scan ',num2str(i)])
    ylabel('Occurrence (n)')
    box off

end
xlabel('Locations (deg)')


%graph prior in Main trials for each scan, second contrast
figure('color','w','position',[1236 841 560 259]);
SLgraphBigTitle(['Prior for Main trials, each scan, contrast ',num2str(uniqCon(2))])

for i = 1 : nScanTot
    
    subplot(nScanTot,1,i)
    
    %Main trials and first contrast condition
    MainTrials = task{i}.parameter.loc.trialType==1;
    thisCon = task{i}.parameter.loc.con==uniqCon(2) & MainTrials==1;
    
    hist(task{i}.parameter.loc.series(thisCon),tmpMain.parameter.loc.sample.degree)
    xlim([5 355])
    title(['Scan ',num2str(i)])
    ylabel('Occurrence (n)')
    box off

end
xlabel('Locations (deg)')


%graph prior in Main trials for each scan, second contrast
figure('color','w','position',[1236 841 560 259]);
SLgraphBigTitle(['Prior for Catch trials, each scan, contrast ',num2str(uniqCon(1))])

for i = 1 : nScanTot
    
    subplot(nScanTot,1,i)
    
    %Main trials and first contrast condition
    CatchTrials = task{i}.parameter.loc.trialType==0;
    thisCon = task{i}.parameter.loc.con==uniqCon(1) & CatchTrials==1;
    
    hist(task{i}.parameter.loc.series(thisCon),tmpMain.parameter.loc.sample.degree)
    xlim([5 355])
    title(['Scan ',num2str(i)])
    ylabel('Occurrence (n)')
    box off

end
xlabel('Locations (deg)')

%graph prior in Main trials for each scan, second contrast
figure('color','w','position',[1236 841 560 259]);
SLgraphBigTitle(['Prior for Catch trials, each scan, contrast ',num2str(uniqCon(2))])

for i = 1 : nScanTot
    
    subplot(nScanTot,1,i)
    
    %Main trials and first contrast condition
    CatchTrials = task{i}.parameter.loc.trialType==0;
    thisCon = task{i}.parameter.loc.con==uniqCon(2) & CatchTrials==1;
    
    hist(task{i}.parameter.loc.series(thisCon),tmpCatch.parameter.loc.sample.degree)
    xlim([5 355])
    title(['Scan ',num2str(i)])
    ylabel('Occurrence (n)')
    box off

end
xlabel('Locations (deg)')


%save task for each scan
allScansTask = task;
clear task

%if task scan exists find last task scan and increment from there,
%otherwise start from 0.
d = dir('taskScan*');
if ~isempty(d)
    lastScan = str2num(d(end).name(length('taskScan')+1:length('taskScan')+2));
    fprintf(['Task parameters were found in the directory. I created new scans starting from taskScan0', num2str(lastScan),'\n'])
else
    lastScan = 0;
end

%save scans counting from last scan
fprintf('----------------------------------------- \n')
fprintf('I am now saving your scans....please wait \n')
fprintf('----------------------------------------- \n')

for i = 1: length(allScansTask)
    task = allScansTask{i};
    
    %clean format
    if lastScan+i>=10
        scanName = 'taskScan';
    else 
        scanName = 'taskScan0';
    end
    
    %format prior
    myPrior = num2str(Pk);
    myPrior(myPrior=='.') = '_';
    scanName = [scanName,num2str(lastScan+i),'_Prior',myPrior];
    save(scanName,'task')
    
    %status
    fprintf([scanName,' was saved. \n'])
end
fprintf('done \n')






