

<<<<<<< HEAD

% slloadClassifAccu.m
%
%load data created by "slfmriClassify2.m"

%move to data path
cd /Volumes/DroboBKUP/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/V1

%load accuracy
load ClassifStckSessmySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06_V1_21-Feb-2016.mat
=======
%slloadClassifAccu.m
%
%
% author: steeve laquitaine
%purpose: load data created and saved by "slfmriClassify2.m"


%rois
rois = {'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'};

%collect decoding accuracies
for i = 1 : length(rois)
    
    roi = rois{i};
    
    %move to data path
    %cd(['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/' roi])
    %cd(['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/mySwitch=1_x_myRandomCoh=0.12mySwitch=2_x_myRandomCoh=0.12/accAtTimeleaveOneOutfisherbalancByRemovI/' roi])
    %cd(['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/' roi])
    %cd(['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06/accAtTimeleaveOneOutfisherbalancByRemovI/' roi])
    cd(['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomDir_x_mySwitch=1/accAtTimeleaveOneOutfisherbalancByRemovI/' roi])  
    
    %load accuracy
    % load(['ClassifStckSessmySwitch=1_x_myRandomCoh=0.06mySwitch=2_x_myRandomCoh=0.06_' roi '_21-Feb-2016.mat'])
    % load(['ClassifStckSessmySwitch=1_x_myRandomCoh=0.12mySwitch=2_x_myRandomCoh=0.12_' roi '_21-Feb-2016.mat'])
    % load(['ClassifStckSessmyRandomDir_x_mySwitch=1_x_myRandomCoh=0.06_' roi '_22-Feb-2016.mat'])
    % load(['ClassifStckSessmyRandomDir_x_mySwitch=2_x_myRandomCoh=0.06_' roi '_22-Feb-2016.mat'])
    load(['ClassifStckSessmyRandomDir_x_mySwitch=1_' roi '_22-Feb-2016.mat'])
    
    acc(i) = round(c.myClasf.raw.fullSg.correct*100)/100;
    ste(i) = round(c.myClasf.raw.fullSg.correctSTE*100)/100;
    %c.myClasf.raw.TheoretChance
      
end
>>>>>>> 608bdce04de4d0b0a6a30ef97e58ba3016b21d1b
