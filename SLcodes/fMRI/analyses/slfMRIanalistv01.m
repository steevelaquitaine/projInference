
%slfMRIanalistv01.m



%% Add Switching variable to data
%151229 (done)
%initialize "switched to prior/dir" variable
%create switching variable for classification
o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';
o.sessPath{1} = [o.rootpath '/s02520150814'];
% o.sessPath{2} = [o.rootpath '/s002520150923'];
% o.sessPath{3} = [o.rootpath '/s002520150925'];
o.myGroup = 'Concatenation';
o.taskNum = 2;
for sess = 1 : length(o.sessPath)
    fprintf('%s%s \n','SESSION:',o.sessPath{1})
    fprintf('%s \n','-----------------------')    
    o.myPath = o.sessPath{sess};
    tmp = slmakeSwitchingVar(o);
    o.behBySess{sess} = tmp;
end

%% Initialize parameters
%Available rois are :
% 'V1','V2','V3','V3A','hV4','MT','V7','IPS','parietalBA39','FEF','vlPFC',
% 'dlPFC','vmPFC','OFC','aPFC'
%make sure ws always clear !
clear
o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';
o.stckPath = o.rootpath;
o.sessPath{1} = [o.rootpath 's02520150814'];
% o.sessPath{2} = [o.rootpath 's002520150923'];
% o.sessPath{3} = [o.rootpath 's002520150925'];
o.myGroup = 'Concatenation';
o.myScan = 1;
o.myROIname = {'V1'};%,'V2','V3','V3A','hV4','MT','V7','IPS','parietalBA39','FEF','vlPFC','dlPFC','vmPFC','OFC','aPFC'};
o.taskNum = 2;
o.phaseNum = 1;
o.segmentNum = 2;   
o.anatFileName = 's0025_flatL_WM_occipital_Rad90';
o.anatFilePath = '~/data/datafMRI/mlrAnatDB/s0025/mlrBaseAnatomies/';
o.myROIpath = '~/data/datafMRI/mlrAnatDB/s0025/mlrROIs/';
% o.erAnalFile = 'erAnalMotionByAll.mat';



%% ----------------------------------- Classification Analyses ---------------------------
%% classify directions
slfmriClassify2(o,'myRandomDir','accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');
cd slStckAnalyses/Concatenation/classif/myRandomDir/accAtTimeleaveOneOutfisherbalancByRemovI/ 

%% classify coh
slfmriClassify2(o,'myRandomCoh','accuracyAtTime',[7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');
cd slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/ 

%% classify switching (done)
myClass = {{'mySwitch=1'},{'mySwitch=2'}};
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');
cd slStckAnalyses/Concatenation/classif/mySwitch=1mySwitch=2/accAtTimeleaveOneOutfisherbalancByRemovI/ 



%% ----------------------------- Classification Analyses with permutation ---------------------------
%between 5 hours and 7 by roi (500 perm)
slfmriClassify2(o,'myRandomDir','accuracyAtTime',[7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher',...
    'balancByRemovI=1','chanceByPerm','nPerm=500','permutationBal=1');
% cd slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal/

%% classify switching
slfmriClassify2(o,{{'mySwitch=1'},{'mySwitch=2'}},'accuracyAtTime',...
    [7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1','chanceByPerm','nPerm=500','permutationBal=1');
cd slStckAnalyses/Concatenation/classif/mySwitch=1mySwitch=2/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal/
%p = length(find(c.myClasf.raw.fullSgShuf.corrects > c.myClasf.raw.fullSg.correct))/length(c.myClasf.raw.fullSgShuf.corrects)*100;

%% classify coh
slfmriClassify2(o,'myRandomCoh','accuracyAtTime',...
    [7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1','chanceByPerm','nPerm=500','permutationBal=1');
cd slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal/


%% ----------------------------------- Dealing with confounds ---------------------------

%% classify switching by coherence (to remove confound)
%we might not decode switching but coherence because subjects switch more
%to prior at low coherence. So we decode switching separaetely by coherence
%% coh 0.06
myClass = {{'mySwitch=1_x_myRandomCoh=0.06'},{'mySwitch=2_x_myRandomCoh=0.06'}};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');
cd slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/

%% coh 0.12
myClass = {{'mySwitch=1_x_myRandomCoh=0.12'},{'mySwitch=2_x_myRandomCoh=0.12'}};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');
cd slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/

%% permutation (leaveOneOut, 1 h?)
%see loadChancePermutationData.m for data summary
myClass = {{'mySwitch=1_x_myRandomCoh=0.06'},{'mySwitch=2_x_myRandomCoh=0.06'}};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher',...
    'balancByRemovI=1','chanceByPerm','nPerm=500','permutationBal=1');
cd slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/

%% permutation (kFold, faster (9 min/roi))
%see loadChancePermutationData.m for data summary
%% coh 0.06
%10 folds x 10 instances (9min)
myClass = {{'mySwitch=1_x_myRandomCoh=0.06'},{'mySwitch=2_x_myRandomCoh=0.06'}};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','kFoldCV','fisher',...
    'balancByRemovI=1','chanceByPerm','nPerm=500','permutationBal=1');
%cd slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimekFoldCVfisherbalancByRemovI/

%% 0.12
myClass = {{'mySwitch=1_x_myRandomCoh=0.12'},{'mySwitch=2_x_myRandomCoh=0.12'}};   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','kFoldCV','fisher',...
    'balancByRemovI=1','chanceByPerm','nPerm=500','permutationBal=1');
%cd slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.12/accAtTimekFoldCVfisherbalancByRemovI/



%% ----------------------------------- Neural mecha of switching ---------------------------

%% classify directions when switch to prior
%Subjects might switch to prior because evidence is flat (Bayes) in which
%case we should not be able to decode motion directions
%coh 0.06
%#number instances for each class
%        sessions  total
%15  --> 3  2  3    8
%85  --> 5  2  2    9
%155 --> 22 8  15  45
%225 --> 0  0  0    0
%295 --> 22 5  15  52
%Dataset was balanced at 8 instances per class
myClass = 'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06';   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

%% coh 0.12
%Cancelled analysis: too few samples
%#number instances for each class
%        sessions  total
%15  --> 1  0  1    2
%85  --> 0  0  0    0
%155 --> 6  0  4   10
%225 --> 0  0  0    0
%295 --> 8  0  5   13
% myClass = 'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.12';   
% slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
%     'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');

%% classify directions when switch to direction
%more of a control to check that oss of accuracy when subjects switched to
%prior is not due to smaller sample size but is a real effect. Here we
%should get good accuracy
%coh 0.06
%#number instances for each class
%        sessions  total
%15  --> 11  10  13   34
%85  --> 11  10  13   34
%155 --> 17  12  20   49
%225 -->  0   0   0    0
%295 --> 17  17  21   55
%Dataset was balanced at 32 instances per class
myClass = 'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06';   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');


%% classify directions when switch to direction (pooled coh)
%Coherence and directions have no reason to be correlated when chosing
%prior so we pool all coherences to get more samples
%
%pooled coherences
%#number instances for each class
%        sessions  total
%15  -->  4   2  4  10
%85  -->  5   2  2   9
%155 -->  28  8  19  55
%225 -->  0   0  0    0
%295 -->  30  5  20  55
%Dataset was balanced at 9 instances per class
myClass = 'myRandomDir_x_mySwitch=1';   
slfmriClassify2(o,myClass,'accuracyAtTime',[7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14;7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');








%-------------------- Voxel population activity analysis ------------------

%% create and add switching variable
%Only need to run once (done)
%initialize "switched to prior [1]/dir[2]" variable
%create and add switching variable
slfmriStckAddSwitchingVar

%% initialize analysis
slfmriInitAnalysisTaskDotDirfMRI05

%% get database (responses and variables) stacked over sessions
%calculate trial-responses to motion
o.myROIname = {'V1'};%{'V1','V2','V3','V3A','hV4','MT','V7','IPS','parietalBA39','FEF','vlPFC','dlPFC','vmPFC','OFC','aPFC'};
slfmriClassify2(o,'myRandomCoh','accuracyAtTime',[7 14],...
    'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1');
cd(['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/' o.myROIname{1}])
%create stacked database of instances and associated variables
[o,sessPath] = slfmriInitGetInstancedb('V1','ClassifStckSessmyRandomCoh_V1.mat');
slfmriStckGetInstancedb

%% check voxels goodness of fit and params
% slfMRIcheckVoxels

%% only use selected voxels
%select voxels for subsequent analyses
%Criteria are: > 1% signal change between and less/more preferred direction
%and high r-squared
%get voxel tuning params and goodness of fit
load d
o.roiname = 'V1';
o.dataPath = ['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/' o.roiname '/'];
o = slfmriGetVoxPopActivity2('myRandomDir',{'mySwitch=1','mySwitch=2'},'myRandomCoh=0.06',85,o,'useinputdb',d,'getVoxTuning','vmfitr2cutoff=0');
load VoxelsParams.mat
%select best voxels (r2 and width threshold)
R2thresh = 0.9;
kThresh = 0.0049;
%best voxels
voxFiltR2_Knum = find(rsquared1>=R2thresh & ks1>=kThresh);
%best widths
ks1FiltR2_K = ks1(rsquared1>=R2thresh & ks1>=kThresh);
dfilt = d;
dfilt.instances = d.instances(:,voxFiltR2_Knum);
d.voxnum = voxFiltR2_Knum;
description = {'voxels were selected based on high R2 and strong selectivity (> 1% signal change between less/more prefdirection'};
save('dfilt','d','description')
fprintf('Strongest tuning r-squared: \n');disp(rsquared1(voxFiltR2_Knum))
fprintf('Strongest tuning widths with 0.9 R2 (vm K): \n');disp(ks1FiltR2_K)


%% plot selected voxel pop of activity by switching
%this gets voxel population (response of voxels sorted by their selectivities)
%sorted by "switching-to-prior" or "switching-to-direction".
%population of activity displays the average response (y-axis) of each voxel 
%(x-axis) to a displayed motion direction. Responses are first averaged
%over repeats within voxels then averaged over voxels with same
%selectivity.  Voxel selectivities are evaluated by fitting voxels with von mises and
%taking the argmax direction.
%The instances are in percent signal change and have been set to have a
%mean of 1 during concatenation

%prior from switch-to-direction to switch-to-prior
%note: speed: parfor (12 cores): 0.14 sec per voxel
%status: tested !!!
o.roiname = 'V1';
o.dataPath = ['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/' o.roiname '/'];
o = slfmriGetVoxPopActivity2('myRandomDir',{'mySwitch=1','mySwitch=2'},'myRandomCoh=0.06',85,o,'useinputdb',dfilt,'getVoxTuning','vmfitr2cutoff=0');
SLmergeFigures([3 4])

%% voxel pop activity by coherence
%this gets voxel population (response of voxels sorted by their selectivities)
%sorted by coherence. Population of activity displays the average response 
%(y-axis) of each voxel (x-axis) to a particular motion direction. Voxel 
%selectivities are evaluated by fitting voxels with von mises and
%taking the argmax direction.
% o.roiname = 'V1';
% o.dataPath = ['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/' o.roiname '/'];
% o = slfmriGetVoxPopActivity2('myRandomDir',{'myRandomCoh=0.06','myRandomCoh=0.12'},'mySwitch=2',85,o,'useinputdb',dfilt,'getVoxTuning','vmfitr2cutoff=0');


%note: we might want to cross-validate selectivities

%% distribution of all voxel selectivities
load VoxelsParams
figure('color','w')
hold all
subplot(121)
hist(modes1,1:10:360)
slsetHistColor([.5 .5 .5],'none')
vline(225,':b')
box off
xlim([0 360]); title('Switch-to-prior')
subplot(122)
hist(modes2,1:10:360)
slsetHistColor([.5 .5 .5],'none')
vline(225,':b')
box off
xlim([0 360]);title('Switch-to-dir')

%% plot population activity by switching together
SLmergeFigures([3 4])

%% Check the tuning fit of individual voxel
%responses are missing at prior mean because it is impossible to classify
%estimates as switching to prior or direction when motion direction is
%displayed at prior mean. Thus those trials are removed.
%status: tested !!!
for voxnum = 26%:299
    st = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch','myRandomCoh=0.12',d,'fminsearch');
    vline(225,'k:')
    pause(0.3)
end













