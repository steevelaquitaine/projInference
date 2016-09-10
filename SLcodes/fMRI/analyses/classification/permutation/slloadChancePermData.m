
%loadChancePermutationData.m
%fisher leaveOneOut on 3 stacked sessions pilot s025
%
%see slPlotClassifAccu.m


%% COHERENCE
%path
%cd ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal
cd IPS

%% data
%load  ClassifStckSessmyRandomCoh_V2_17-Feb-2016.mat
% load ClassifStckSessmyRandomCoh_V3_18-Feb-2016.mat
%load ClassifStckSessmyRandomCoh_V3A_18-Feb-2016.mat
%load ClassifStckSessmyRandomCoh_MT_18-Feb-2016.mat
%load ClassifStckSessmyRandomCoh_hV4_18-Feb-2016.mat
%load ClassifStckSessmyRandomCoh_V7_18-Feb-2016.mat
load ClassifStckSessmyRandomCoh_IPS_19-Feb-2016.mat

%% chance accuracy (permuted labels)
hist(c.myClasf.raw.fullSgShuf.corrects)
xlim([0 1])

%p-value
nPermAboveAccuracy = sum(c.myClasf.raw.fullSgShuf.corrects > c.myClasf.raw.fullSg.correct);
nPerm = length(c.myClasf.raw.fullSgShuf.corrects);
ptmp = nPermAboveAccuracy/nPerm

%% stored p-values
roi = {'V1','V2','V3','V3A','MT','hV4','V7','IPS'};
p = [0 0 0 0 0 0 0 0];




%% Direction
%path
cd ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomDir/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal
cd IPS

%% data
%load ClassifStckSessmyRandomDir_V1_30-Jan-2016.mat
%load ClassifStckSessmyRandomDir_V2_04-Feb-2016.mat
%load ClassifStckSessmyRandomDir_V3_05-Feb-2016.mat
%load ClassifStckSessmyRandomDir_V3A_05-Feb-2016.mat
%load ClassifStckSessmyRandomDir_MT_07-Feb-2016.mat
%load ClassifStckSessmyRandomDir_hV4_05-Feb-2016.mat
%load ClassifStckSessmyRandomDir_V7_08-Feb-2016.mat
load ClassifStckSessmyRandomDir_IPS_08-Feb-2016.mat

%% chance accuracy (permuted labels)
hist(c.myClasf.raw.fullSgShuf.corrects)
xlim([0 1])

%p-value
nPermAboveAccuracy = sum(c.myClasf.raw.fullSgShuf.corrects > c.myClasf.raw.fullSg.correct);
nPerm = length(c.myClasf.raw.fullSgShuf.corrects);
ptmp = nPermAboveAccuracy/nPerm


%% stored p-values
roi = {'V1','V2','V3','V3A','MT','hV4','V7','IPS'};
p =   [0 0 0 0.004 0.22 0 0 0];





%% Switching
%path
%cd ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1mySwitch=2/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal
cd ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI_chanceByPermpermutationBal/
cd V1
% load ClassifStckSessmySwitch=1mySwitch=2_V1_06-Jan-2016.mat
% load ClassifStckSessmySwitch=1mySwitch=2_V2_08-Jan-2016.mat
% load ClassifStckSessmySwitch=1mySwitch=2_V3_08-Jan-2016.mat
% load ClassifStckSessmySwitch=1mySwitch=2_V3A_06-Jan-2016.mat
% load ClassifStckSessmySwitch=1mySwitch=2_MT_07-Jan-2016.mat
% load ClassifStckSessmySwitch=1mySwitch=2_hV4_09-Jan-2016.mat
% load ClassifStckSessmySwitch=1mySwitch=2_V7_09-Jan-2016.mat
%load ClassifStckSessmySwitch=1mySwitch=2_IPS_07-Jan-2016.mat
load ClassifStckSessmySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06_V1.mat
hist(c.myClasf.raw.fullSgShuf.corrects)
xlim([0 1])
%p-value
nPermAboveAccuracy = sum(c.myClasf.raw.fullSgShuf.corrects > c.myClasf.raw.fullSg.correct);
nPerm = length(c.myClasf.raw.fullSgShuf.corrects);
ptmp = nPermAboveAccuracy/nPerm
%stored p-values
roi = {'V1','V2','V3','V3A','MT','hV4','V7','IPS'};
p = [0 ];


%% Switching by coh (dealing with coherence confound)
%10-fold cross validation
%0.06
cd ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimekFoldCVfisherbalancByRemovI_chanceByPermpermutationBal/
roi = {'V1','V2','V3','V3A','MT','hV4','V7','IPS'};
for i = 1 : length(roi)
    cd(roi{i})
    load(['ClassifStckSessmySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06_' roi{i} '.mat'])
    fprintf('percent correct: \n'); disp(c.myClasf.raw.fullSg.correct)
    %hist(c.myClasf.raw.fullSgShuf.corrects)
    %xlim([0 1])
    %p-value
    nPermAboveAccuracy = sum(c.myClasf.raw.fullSgShuf.corrects > c.myClasf.raw.fullSg.correct);
    nPerm = length(c.myClasf.raw.fullSgShuf.corrects);
    ptmp = nPermAboveAccuracy/nPerm;
    %stored
    p(i) = ptmp; 
    acc(i) = round(c.myClasf.raw.fullSg.correct*100)/100*100;
    ste(i) = round(c.myClasf.raw.fullSg.correctSTE*100)/100*100;
    cd ..
end


%% 0.12
cd ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.12mySwitch=2_x_myRandomCoh=0.12/accAtTimekFoldCVfisherbalancByRemovI_chanceByPermpermutationBal/
roi = {'V1','V2','V3','V3A','MT','hV4','V7','IPS'};
for i = 1 : length(roi)
    cd(roi{i})
    load(['ClassifStckSessmySwitch=1_x_myRandomCoh=0.12mySwitch=2_x_myRandomCoh=0.12_' roi{i} '.mat'])
    fprintf('percent correct: \n'); disp(c.myClasf.raw.fullSg.correct)
    %hist(c.myClasf.raw.fullSgShuf.corrects)
    %xlim([0 1])
    %p-value
    nPermAboveAccuracy = sum(c.myClasf.raw.fullSgShuf.corrects > c.myClasf.raw.fullSg.correct);
    nPerm = length(c.myClasf.raw.fullSgShuf.corrects);
    ptmp = nPermAboveAccuracy/nPerm;
    %stored
    p(i) = ptmp; 
    acc(i) = round(c.myClasf.raw.fullSg.correct*100)/100*100;
    ste(i) = round(c.myClasf.raw.fullSg.correctSTE*100)/100*100;
    cd ..
end



