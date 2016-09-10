
%SLinitTrainingTask
%
%     author: steeve laquitaine
%       date: 140505
%    purpose: Create a behavioral training run by concatenating fMRI scans parameters
%
%      usage:
% 
%         %80 deg prior and 40 deg prior
%         task = SLinitTrainingTask('training_Prior0_74848',[.5 0 0],...
%             {'taskScan01_Prior0_74848',...
%             'taskScan02_Prior0_74848',...
%             'taskScan03_Prior0_74848',...
%             'taskScan04_Prior0_74848',...
%             'taskScan05_Prior0_74848'})
% 
%         task = SLinitTrainingTask('training_Prior2_7714',[1 0.6 0],...
%             {'taskScan06_Prior2_7714',...
%             'taskScan07_Prior2_7714',...
%             'taskScan08_Prior2_7714',...
%             'taskScan09_Prior2_7714',...
%             'taskScan10_Prior2_7714'})

%Description: this code uses the mat files created by SLinitEachScanParams
%             to generate training .mat files. 
%             To make parameter changes, first run SLinitEachScanParams
%             then run SLinitTrainingTask. SLinitTrainingTask will inherit
%             the parameters of SLinitEachScanParams.


function task = SLinitTrainingTask(nmParams,colors,SubRunParams)

%move to TaskParameters directory
myPath = mfilename('fullpath');
cd(fileparts(myPath))
cd ..

nSubRuns = length(SubRunParams);

%initialize parameters
           task.parameter.loc.series = [];
              task.parameter.loc.con = []; 
    task.parameter.loc.uniqtrialType = []; 
         task.parameter.loc.strength = []; 
             task.parameter.loc.mean  = []; 
           task.parameter.loc.sample = []; 
          task.parameter.loc.sampsiz = []; 
           task.parameter.loc.ntrial = []; 

%concatenate SubRuns params
for i = 1 : nSubRuns
    
    tmp{i} = load(['TaskParameters/',SubRunParams{i},'.mat']);    
    
           task.parameter.loc.series = [task.parameter.loc.series        tmp{i}.task.parameter.loc.series ]; 
              task.parameter.loc.con = [task.parameter.loc.con           tmp{i}.task.parameter.loc.con]; 
         task.parameter.loc.strength = unique([task.parameter.loc.strength      tmp{i}.task.parameter.loc.strength]); 
             task.parameter.loc.mean = unique([task.parameter.loc.mean   tmp{i}.task.parameter.loc.mean]); 
           task.parameter.loc.sample = unique([task.parameter.loc.sample tmp{i}.task.parameter.loc.sample]); 
           task.parameter.loc.ntrial = [task.parameter.loc.ntrial        tmp{i}.task.parameter.loc.ntrial]; 
    
end

task.parameter.loc.ntrial = length(task.parameter.loc.series);
task.parameter.loc.sampsiz = length(task.parameter.loc.sample);
task.parameter.loc.sampleCon = unique(task.parameter.loc.con);
task.parameter.loc.nCon = length(task.parameter.loc.sampleCon );
task.parameter.loc.trialType = zeros(task.parameter.loc.ntrial,1);
task.parameter.loc.uniqtrialType = 0;

for i = 1 : task.parameter.loc.nCon
    task.parameter.loc.TrueTrialnumperCon(i) = sum(task.parameter.loc.con == task.parameter.loc.sampleCon(i));
end

%save
clearex('params')

%move to TaskParameters directory
myPath = mfilename('fullpath');
cd(fileparts(myPath))
cd ../TaskParameters



%graph prior

%polar and hist distributions
%Draw polar distributions
figure('color','w')
myOdd = SLmakeOdd(0,task.parameter.loc.nCon*2);

for j = 1 : task.parameter.loc.nCon
    subplot(task.parameter.loc.nCon,2,myOdd(j));
    
    for i = 1 : task.parameter.loc.sampsiz
        task.parameter.loc.countFeat_and_StimStrength(j,i) = sum(task.parameter.loc.series==task.parameter.loc.sample(i) & task.parameter.loc.con==task.parameter.loc.sampleCon(j));
    end
    SLcircHist(task.parameter.loc.sample,task.parameter.loc.countFeat_and_StimStrength(j,:),colors,5)
end

%Draw discrete distributions (linear space)
fig.name = nmParams;
myEven = SLmakeEven(0,task.parameter.loc.nCon*2);
for i = 1 : task.parameter.loc.nCon
    
    subplot(task.parameter.loc.nCon,2,myEven(i));
    
    %experimental distribution
    SLdrawhistdir(task.parameter.loc.sample,task.parameter.loc.countFeat_and_StimStrength(i,:),fig,colors);
    
    %overlap continuous von Mises
    [~,count] = SLmakeDicreteCountVMdensity(task.parameter.loc.sample,task.parameter.loc.mean,task.parameter.loc.strength,task.parameter.loc.TrueTrialnumperCon(i));
    plot(task.parameter.loc.sample,count,'color','k');
    text(0.9*max(task.parameter.loc.sample),.90*max(task.parameter.loc.countFeat_and_StimStrength(i,:)),['con ' num2str(task.parameter.loc.sampleCon (i))],'fontsize',14)
    
    %info
    title(['Prior k : ' num2str(task.parameter.loc.strength),', trial nb : ' num2str(task.parameter.loc.TrueTrialnumperCon(i))])
    lg = legend('Actual disc. dis','','Orig.cont.dis','location','NorthWest');
    set(lg,'fontsize',10)
    legend('boxoff')
    
    %graphics
    task.parameter.loc.contDistInCount(i,:) = count;
    maxp = max([task.parameter.loc.countFeat_and_StimStrength(i,:) task.parameter.loc.contDistInCount(i,:)]);
    ylim([0 maxp])
    
end
xlabel({'Displayed Feature (deg)'})



SLautobackup(gcf,nmParams,'.fig')
SLautobackup(task,nmParams,'.mat')