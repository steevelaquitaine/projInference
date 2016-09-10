
%slfMRIclassifSAdisplay.m
%
%
% author: steeve laquitaine
%   date: 151016
%purpose: plot summary results from classification of coherence from
%         MT,V3A,IPS and V1 produced by slfMRIclassifAnalSensoryAreas2.m
%  usage:
%
%
%         o.rootpath = '/Volumes/Transcend/data/sltaskdotdirfmri05/'; %classif data path
%         o.myAnalysis = {'accAtTime','kFoldCV','_chanceByPerm'};
%         o.myClass = 'mySwitch=1_x_myRandomDir=15';
%         o.nRois = 2;
%         o.stk   = 1; %aligned sessions
%
%         [mx,mxz] = slfMRIclassifSAdisplay(o)


function [mx,mxz] = slfMRIclassifSAdisplay(o)

%check cell
if ~iscell(o.myAnalysis)
    o.myAnalysis = {o.myAnalysis};
end

figure('color','w')

%myRandomCoh mixte dir. leaveOneOut
if all(strcmp(o.myClass,'myRandomCoh')) && all(strcmp(o.myAnalysis,'accAtTime'))
    [mx,mxz] = randCohMixteDir;
end
%myRandomCoh mixte dir r2cutoff leaveOneOut
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimer2cutoff')
    [mx,mxz] = myRandomCohaccAtTimer2cutoff;
end
%myRandomCoh mixte dir leaveOneKmeaninstanceout
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimeslleaveOneOfkmeanInstanceOut')
    [mx,mxz] = myRandomCohaccAtTimeslleaveOneOfkmeanInstanceOut;
end
%myRandomCoh mixte dir r2cutoff leaveOneKmeaninstanceout
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimer2cutoffslleaveOneOfkmeanInstanceOut')
    [mx,mxz] = myRandomCohaccAtTimer2cutoffslleaveOneOfkmeanInstanceOut;
end
%myRandomCoh for each dir leaveOneOut
if length(o.myClass)==26 && strcmp(o.myClass(1:26),'myRandomCoh_x_myRandomDir=') && all(strcmp(o.myAnalysis,'accAtTime'))
    [mx,mxz] = myRandomCohEachDiraccAtTimeleaveOneOut;
end
%myRandomCoh for each dir r2cutoff leaveOneOut
if  length(o.myClass)==26 && strcmp(o.myClass(1:26),'myRandomCoh_x_myRandomDir=') && strcmp([o.myAnalysis{:}],'accAtTimer2cutoff')
    [mx,mxz] = myRandomCohEachDiraccAtTimer2cutoffleaveOneOut;
end
%myRandomCoh for each dir leaveOneKmeaninstanceout
if length(o.myClass)==26 && strcmp(o.myClass(1:26),'myRandomCoh_x_myRandomDir=') && strcmp([o.myAnalysis{:}],'accAtTimeslleaveOneOfkmeanInstanceOut')
    [mx,mxz] = myRandomCohEachDiraccAtTimeslleaveOneOfkmeanInstanceOut;
end
%myRandomCoh for each dir r2 cutoff leaveOneKmeaninstanceout
if  length(o.myClass)==26 && strcmp(o.myClass(1:26),'myRandomCoh_x_myRandomDir=') && strcmp([o.myAnalysis{:}],'accAtTimer2cutoffslleaveOneOfkmeanInstanceOut')
    [mx,mxz] = myRandomCohEachDiraccAtTimesr2cutofflleaveOneOfkmeanInstanceOut;
end
%myRandomCoh for regressed-out-direction-weights leaveOneOut
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimeleaveOneOutregOutVar2')
    [mx,mxz] = myRandomCohaccAtTimeleaveOneOutregOutVar2;
end
%myRandomCoh for regressed-out-direction-weights r2cutoff leaveOneOut
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimer2cutoffleaveOneOutregOutVar2')
    [mx,mxz] = accAtTimer2cutoffleaveOneOutregOutVar2;
end
%myRandomCoh for regressed-out-direction-weights leave one of k mean instances out
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimeslleaveOneOfkmeanInstanceOutregOutVar2')
    [mx,mxz] = accAtTimeslleaveOneOfkmeanInstanceOutregOutVar2;
end
%myRandomCoh for regressed-out-direction-weights r2cutoff leave one of k mean instances out
if all(strcmp(o.myClass,'myRandomCoh')) && strcmp([o.myAnalysis{:}],'accAtTimer2cutoffslleaveOneOfkmeanInstanceOutregOutVar2')
    [mx,mxz] = accAtTimer2cutoffslleaveOneOfkmeanInstanceOutregOutVar2;
end
%myRandomCoh for each dir leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:26),'myRandomCoh_x_myRandomDir=') && all(strcmp(o.myAnalysis,'accAtTime')) && o.stk==1
    [mx,mxz] = myRandomCohEachDiraccAtTimeleaveOneOutStk;
end
%myRandomDir for each coh leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:26),'myRandomDir_x_myRandomCoh=') && all(strcmp(o.myAnalysis,'accAtTime')) && o.stk==1
    [mx,mxz] = myRandomDirEachCohaccAtTimeleaveOneOutStk(o);
end
%mySwitch for each coh leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:25),'mySwitch=1_x_myRandomCoh=') && all(strcmp(o.myAnalysis,'accAtTime')) && o.stk==1
    [mx,mxz] = mySwitchbyCohaccAtTimeleaveOneOutStk(o);
end
%mySwitch for each coh leaveOneOut (stacked sess)
if all(strcmp(o.myClass,{'mySwitch=1','mySwitch=2'})) && all(strcmp(o.myAnalysis,'accAtTime')) && o.stk==1
    [mx,mxz] = mySwitchaccAtTimeleaveOneOutStk(o);
end
%mySwitch for each coh leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:25),'mySwitch=1_x_myRandomDir=') && all(strcmp([o.myAnalysis{:}],'accAtTimekFoldCV')) && o.stk==1
    [mx,mxz] = mySwitchbyDiraccAtTimeleaveOneOutStk(o);
end
%mySwitch for each coh leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:25),'mySwitch=1_x_myRandomDir=') && all(strcmp([o.myAnalysis{:}],'accAtTimekFoldCV_chanceByPerm')) && o.stk==1
    [mx,mxz] = mySwitchbyDirWithPermaccAtTimeleaveOneOutStk(o);
end
%mySwitch for each coh leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:25),'mySwitch=1_x_myRandomCoh=') && all(strcmp([o.myAnalysis{:}],'accAtTimekFoldCVsvm')) && o.stk==1
    [mx,mxz] = mySwitchbyCohaccAtTimeleaveOneOutStk(o);
end
%mySwitch for each coh leaveOneOut (stacked sess)
if ~iscell(o.myClass) && strcmp(o.myClass(1:25),'mySwitch=1_x_myRandomCoh=') && all(strcmp([o.myAnalysis{:}],'accAtTimekFoldCVsvm_chanceByPerm')) && o.stk==1
    [mx,mxz] = mySwitchbyCohWithPermaccAtTimeleaveOneOutStk(o);
end

function [mx,mxz] = randCohMixteDir

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTime/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh mixted dir Attime')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
function [mx,mxz] = myRandomCohaccAtTimer2cutoff

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimer2cutoff/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh mixted dir r2cutoff Attime')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
xlabel('Raw/Z-scored over vols')
function [mx,mxz] = myRandomCohaccAtTimeslleaveOneOfkmeanInstanceOut

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimeslleaveOneOfkmeanInstanceOut/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh mixted dir Attime leaveOneOfkmeanInstanceOut')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
xlabel('Raw/Z-scored over vols')
function [mx,mxz] = myRandomCohaccAtTimer2cutoffslleaveOneOfkmeanInstanceOut

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimer2cutoffslleaveOneOfkmeanInstanceOut/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh mixted dir Attime leaveOneOfkmeanInstanceOut')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
function [mx,mxz] = myRandomCohEachDiraccAtTimeleaveOneOut

clear
rois = {'MT','V3A','IPS','V1'};
nRois = length(rois);
di  = [15 85 155 225 295];
nDir = length(di);

%cd to data saved for each direction
%collect accuracies and stds.
for d = 1 : nDir
    analpath{d} = sprintf('%s%i%s','~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh_x_myRandomDir=',di(d),'/accAtTime/');
    for i = 1 : nRois
        
        %cd path
        cd([analpath{d} rois{i}]);
        
        %get file
        myfile = dir('*mat');
        clear('c')
        load(myfile.name)
        
        %store accuracy and STE
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;    %raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;    %zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
    end
end


%calculate means of means and ste of means over dirs
%for each roi
%raw
[cbyRoi,cbyRoiste] = slMakeMeanAndStd(mx,mxste);

%z-scored
[cbyRoiz,cbyRoizte] = slMakeMeanAndStd(mxz,mxzste);

%bar
cM    = [cbyRoi cbyRoiz];
csteM = [cbyRoiste cbyRoizte];
xloc  = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('Mean acc. for coh over each dir Attime leaveOneOut')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM' csteM'],[rois rois],{'%c','ste'})
SLprintTable([cM' csteM'],[rois rois],{'%c','ste'})
function [mx,mxz] = myRandomCohEachDiraccAtTimer2cutoffleaveOneOut

clear
rois = {'MT','V3A','IPS','V1'};
nRois = length(rois);
di  = [15 85 155 225 295];
nDir = length(di);

%cd to data saved for each direction
%collect accuracies and stds.
for d = 1 : nDir
    analpath{d} = sprintf('%s%i%s','~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh_x_myRandomDir=',di(d),'/accAtTimer2cutoff/');
    for i = 1 : nRois
        
        %cd path
        cd([analpath{d} rois{i}]);
        
        %get file
        myfile = dir('*mat');
        clear('c')
        load(myfile.name)
        
        %store accuracy and STE
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;    %raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;    %zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
    end
end


%calculate means of means and ste of means over dirs
%for each roi
%raw
[cbyRoi,cbyRoiste] = slMakeMeanAndStd(mx,mxste);

%z-scored
[cbyRoiz,cbyRoizte] = slMakeMeanAndStd(mxz,mxzste);

%bar
cM    = [cbyRoi cbyRoiz];
csteM = [cbyRoiste cbyRoizte];
xloc  = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('Mean acc. for coh over each dir Attime r2cutoff leaveOneOut')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM' csteM'],[rois rois],{'%c','ste'})
function [mx,mxz] = myRandomCohEachDiraccAtTimeslleaveOneOfkmeanInstanceOut

clear
rois = {'MT','V3A','IPS','V1'};
nRois = length(rois);
di  = [15 85 155 225 295];
nDir = length(di);

%cd to data saved for each direction
%collect accuracies and stds.
for d = 1 : nDir
    analpath{d} = sprintf('%s%i%s','~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh_x_myRandomDir=',di(d),'/accAtTimeslleaveOneOfkmeanInstanceOut/');
    for i = 1 : nRois
        
        %cd path
        cd([analpath{d} rois{i}]);
        
        %get file
        myfile = dir('*mat');
        clear('c')
        load(myfile.name)
        
        %store accuracy and STE
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;    %raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;    %zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
    end
end


%calculate means of means and ste of means over dirs
%for each roi
%raw
[cbyRoi,cbyRoiste] = slMakeMeanAndStd(mx,mxste);

%z-scored
[cbyRoiz,cbyRoizte] = slMakeMeanAndStd(mxz,mxzste);

%bar
cM    = [cbyRoi cbyRoiz];
csteM = [cbyRoiste cbyRoizte];
xloc  = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('Mean acc. for coh over each dir Attime leaveOneOfKmeanInstances out')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM' csteM'],[rois rois],{'%c','ste'})
function [mx,mxz] = myRandomCohEachDiraccAtTimesr2cutofflleaveOneOfkmeanInstanceOut

clear
rois = {'MT','V3A','IPS','V1'};
nRois = length(rois);
di  = [15 85 155 225 295];
nDir = length(di);

%cd to data saved for each direction
%collect accuracies and stds.
for d = 1 : nDir
    analpath{d} = sprintf('%s%i%s','~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh_x_myRandomDir=',di(d),'/accAtTimer2cutoffslleaveOneOfkmeanInstanceOut/');
    for i = 1 : nRois
        
        %cd path
        cd([analpath{d} rois{i}]);
        
        %get file
        myfile = dir('*mat');
        clear('c')
        load(myfile.name)
        
        %store accuracy and STE
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;    %raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;    %zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
    end
end


%calculate means of means and ste of means over dirs
%for each roi
%raw
[cbyRoi,cbyRoiste] = slMakeMeanAndStd(mx,mxste);

%z-scored
[cbyRoiz,cbyRoizte] = slMakeMeanAndStd(mxz,mxzste);

%bar
cM    = [cbyRoi cbyRoiz];
csteM = [cbyRoiste cbyRoizte];
xloc  = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('Mean acc. for coh over each dir Attime r2cutoff leaveOneOfKmeanInstances out')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM' csteM'],[rois rois],{'%c','ste'})
function [mx,mxz] = myRandomCohaccAtTimeleaveOneOutregOutVar2

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutregOutVar2/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh regressed-out dir Attime leaveOneOut')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
function [mx,mxz] = accAtTimer2cutoffleaveOneOutregOutVar2

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimer2cutoffleaveOneOutregOutVar2/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh regressed-out dir r2cutoff Attime leaveOneOut')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
xlabel('Raw/Z-scored over vols')
function [mx,mxz] = accAtTimeslleaveOneOfkmeanInstanceOutregOutVar2

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimeslleaveOneOfkmeanInstanceOutregOutVar2/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh regressed-out dir Attime leaveOneofKmeanInstancesOut')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
xlabel('Raw/Z-scored over vols')
function [mx,mxz] = accAtTimer2cutoffslleaveOneOfkmeanInstanceOutregOutVar2

clear
nRois = 4;
rois = {'MT','V3A','IPS','V1'};
analpath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh/accAtTimer2cutoffslleaveOneOfkmeanInstanceOutregOutVar2/';
for i = 1 : nRois
    %cd path
    cd([analpath rois{i}]);
    
    %get file
    myfile = dir('*mat');
    load(myfile.name)
    
    %store accuracy and STE
    mx(i,1)  = c.myClasf.raw.fullSg.correct;
    mx(i,2)  = c.myClasf.raw.fullSg.correctSTE;
    mxz(i,1) = c.myClasf.Zsc.fullSg.correct;
    mxz(i,2) = c.myClasf.Zsc.fullSg.correctSTE;
end

%bar
cbyRoi     = mx(:,1);
cbyRoiste  = mx(:,2);
cbyRoiz    = mxz(:,1);
cbyRoizte  = mxz(:,2);
cM         = [cbyRoi;cbyRoiz];
csteM      = [cbyRoiste;cbyRoizte];
xloc = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('coh regressed-out dir Attime r2cutoff leaveOneofKmeanInstancesOut')

%table
SLprintTable([cM csteM],[rois rois],{'%c','ste'})
xlabel('Raw/Z-scored over vols')
function [mx,mxz] = myRandomCohEachDiraccAtTimeleaveOneOutStk

clear
rois = {'MT','V3A','IPS','V1'};
nRois = length(rois);
di  = [15 85 155 225 295];
nDir = length(di);

%cd to data saved for each direction
%collect accuracies and stds.
for d = 1 : nDir
    %analpath{d} = sprintf('%s%i%s','~/data/datafMRI/sltaskdotdirfmri05/s02520150814/slAnalyses/Concatenation/classif/myRandomCoh_x_myRandomDir=',di(d),'/accAtTime/');
    analpath{d} = sprintf('%s%i%s','~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh_x_myRandomDir=',di(d),'/accAtTime/');
    for i = 1 : nRois
        
        %cd path
        cd([analpath{d} rois{i}]);
        
        %get file
        myfile = dir('*mat');
        clear('c')
        load(myfile.name)
        
        %store accuracy and STE
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;    %raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;    %zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
        
        %get accuracies for shuffled datasets
        for j = 1 : length(c.myClasf.raw.fullSgShuf)
            tp1(j)   = c.myClasf.raw.fullSgShuf(j).correct;    %raw
            tp2(j)   = c.myClasf.Zsc.fullSgShuf(j).correct;    %zscored
        end
        mxS(d,i) = mean(tp1);
        mxzS(d,i) = mean(tp2);
        mxsteS(d,i) = std(tp1);
        mxzsteS(d,i) = std(tp2);
    end
end

%calculate means of means and ste of means over dirs
%for each roi
%raw

% [cbyRoi,cbyRoiste] = slMakeMeanAndStd(mx,mxste);
% [cbyRoiS,cbyRoisteS] = slMakeMeanAndStd(mxS,mxsteS);
%
% %z-scored
% [cbyRoiz,cbyRoizte] = slMakeMeanAndStd(mxz,mxzste);
% [cbyRoizS,cbyRoizteS] = slMakeMeanAndStd(mxzS,mxzsteS);
%calculate means of means and ste of means over dirs
%for each roi
%raw
cbyRoi = mean(mx,1);   %mean over coh
cbyRoiste = std(mx,1); %std over coh
%[cbyRoi,cbyRoiste] = slMakeMeanAndStd(mx,mxste);
cbyRoiS = mean(mxS,1);
cbyRoisteS = std(mxS,1);
%[cbyRoiS,cbyRoisteS] = slMakeMeanAndStd(mxS,mxsteS);

%z-scored
cbyRoiz = mean(mxz,1);
cbyRoizte = std(mxz,1);
%[cbyRoiz,cbyRoizte] = slMakeMeanAndStd(mxz,mxzste);
cbyRoizS = mean(mxzS,1);
cbyRoizteS = std(mxzS,1);
%[cbyRoizS,cbyRoizteS] = slMakeMeanAndStd(mxzS,mxzsteS);

%bar
cM    = [cbyRoi cbyRoiz];
csteM = [cbyRoiste cbyRoizte];
cMS    = [cbyRoiS cbyRoizS]; %shuffled
csteMS = [cbyRoisteS cbyRoizteS];
xloc  = [1:nRois nRois+2:2*nRois+1];
SLdrawBar(cM,[rois rois],xloc,'yError',csteM,'FaceColor',[.8 .8 .8])
hold on; errorbar(xloc,cMS,csteMS)
hline(0.5);
ylabel('Classification accuracy (% correct)')
title('Mean acc. for coh over each dir Attime leaveOneOut')
xlabel('Raw/Z-scored over vols')

%table
SLprintTable([cM' csteM'],[rois rois],{'%c','ste'})
SLprintTable([cMS' csteMS'],[rois rois],{'%cS','steS'})
function [mx,mxz] = mySwitchbyCohWithPermaccAtTimeleaveOneOutStk(o)

rois = o.myROIname;
coi  = [0.06 0.12];
nCoh = length(coi);

%each coherence
fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','Test p-value that %correct comes from chance (null Hyp)')
fprintf('%s \n','-------------------------------------------------------')
for d = 1 : nCoh
    %class
    myClass = ['mySwitch=1_x_myRandomCoh=' num2str(coi(d)) 'mySwitch=2_x_myRandomCoh=' num2str(coi(d))];
    %path
    analpath{d} = [o.rootpath 'slStckAnalyses/' ...
        o.myGroup '/classif/' myClass '/' [o.myAnalysis{:}] '/'];
    
    fprintf('%s %.2f \n','Coh:',coi(d))
    fprintf('%s \n','-------------')
    %each roi
    for i = 1 : o.nRois
        
        %check path
        if exist(analpath{d})
            %data path
            cd([analpath{d} rois{i}]);
        else
            slPrintfStr('slfMRIclassifSAdisplay',' Path does not exist.')
        end
        
        %get file
        myfile = dir('ClassifStckSess*mat');
        clear('c')
        load(myfile.name)
        
        %store %correct and STE
        %actual data
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;%raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;%zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
        
        %get mean %correct for Nx-permuted dataset (chance)
        mxS(d,i) = c.myClasf.raw.fullSgShuf.correct; %raw
        mxzS(d,i) = c.myClasf.Zsc.fullSgShuf.correct; %zscored
        mxsteS(d,i) = c.myClasf.raw.fullSgShuf.correctSTE;
        mxzsteS(d,i) = c.myClasf.Zsc.fullSgShuf.correctSTE;
        
        %95% Confidence interval
        mxCI95S(d,i) = c.myClasf.raw.fullSgShuf.CI95updown(1);
        mxzCI95S(d,i) = c.myClasf.Zsc.fullSgShuf.CI95updown(1);
        
        %calculate p-value
        %statistical significance
        %of each % correct
        %pval is the percent of permutation accuracies
        %above data classification accuracy
        fprintf('%s \n','roi:',o.myROIname{i})
        fprintf('%s \n','---')
        
        %raw
        correctPerm = c.myClasf.raw.fullSgShuf.corrects;%raw
        pval.pval(d,i) = sum(correctPerm > mx(d,i))/length(correctPerm);
        fprintf('%s %.3f \n','p-value (raw):',pval.pval(d,i))
        
        %z-sc
        correctPermz = c.myClasf.Zsc.fullSgShuf.corrects;%z-sc
        pval.pvalz(d,i) = sum(correctPermz > mxz(d,i))/length(correctPermz);
        fprintf('%s %.3f \n','p-value (z-sc):',pval.pvalz(d,i))
        
        %notes
        pval.description = ['p-value for null hypothesis that classification',...
            '%correct comes from chance (N permutations)'];
        save(['pval' o.myROIname{i}],'pval')
    end
end

%----------------- %correct by coh and roi ---------------

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (red (95%CI))'])
    for co = 1 : nCoh
        hold all
        %chance
        bar(2*co-1,mxS(co,roi),'facecolor',[.5 .5 .5])
        %data
        bar(2*co,mx(co,roi),'facecolor',[0 0 0])
        %chance std and Confidence intervals
        %data classif acc set to 0
        %95%CI
        margErr = mxS(co,roi) - mxCI95S(co,roi);
        SLerrorbar([2*co-1 2*co],[mxS(co,roi) mx(co,roi)],...
            'yError',[margErr 0],'MarkerSize=1',...
            'Linestyle=none','Color=r');
    end
    set(gca,'xtick',[1.5 3.5],'xticklabel',coi)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Coherence')
    hline(0.5);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')

%save figure
saveas(gcf, 'figAccData', 'fig')

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (z-scored,red (95%CI))'])
    for co = 1 : nCoh
        hold all
        %chance
        bar(2*co-1,mxzS(co,roi),'facecolor',[.5 .5 .5])
        %data
        bar(2*co,mxz(co,roi),'facecolor',[0 0 0])
        %chance std and Confidence intervals
        %data classif acc set to 0
        %95%CI
        margErrz = mxzS(co,roi) - mxzCI95S(co,roi);
        SLerrorbar([2*co-1 2*co],[mxzS(co,roi) mxz(co,roi)],...
            'yError',[margErrz 0],'MarkerSize=1',...
            'Linestyle=none','Color=r');
        %std
        %SLerrorbar([2*co-1 2*co],[mxzS(co,roi) mxz(co,roi)],...
        %    'yError',[mxzsteS(co,roi) 0],'MarkerSize=1','Linestyle=none')
    end
    set(gca,'xtick',[1.5 3.5],'xticklabel',coi)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Coherence')
    hline(0.5);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')


%save figure
saveas(gcf, 'figAccZscored', 'fig')

%table
for d = 1 : nCoh
    % %correct and ste by roi
    SLprintTable([mx(d,:)' mxste(d,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figAccurAndSTEbyRoi' num2str(d)], 'eps')
    % scored %correct and ste by roi
    SLprintTable([mxz(d,:)' mxzste(d,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figZscAccurAndSTEbyRoi' num2str(d)], 'eps')
    % pvals by roi
    SLprintTable([pval.pval(d,:)' pval.pvalz(d,:)'],rois,{'%c','ste'})
    saveas(gcf, ['pvals' num2str(d)], 'eps')
end
function [mx,mxz] = myRandomDirEachCohaccAtTimeleaveOneOutStk(o)

rois = o.myROIname;
coi  = [0.06 0.12];
nCoh = length(coi);

%each coherence

fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','Test p-value that %correct comes from chance (null Hyp)')
fprintf('%s \n','-------------------------------------------------------')
for d = 1 : nCoh
    
    %analysis path
    analpath{d} = sprintf('%s%s%s',[o.rootpath 'slStckAnalyses/' ...
        o.myGroup '/classif/myRandomDir_x_myRandomCoh=',num2str(coi(d)),'/' [o.myAnalysis{:}] '/']);
    
    fprintf('%s %.2f \n','Coh:',coi(d))
    fprintf('%s \n','-------------')
    %each roi
    for i = 1 : o.nRois
        
        %check path
        if exist(analpath{d})
            %data path
            cd([analpath{d} rois{i}]);
        else
            slPrintfStr('slfMRIclassifSAdisplay',' Path does not exist.')
        end
        
        %get file
        myfile = dir('ClassifStckSess*mat');
        clear('c')
        load(myfile.name)
        
        %store %correct and STE
        %actual data
        %Ndi by Nrois
        mx(d,i)    = c.myClasf.raw.fullSg.correct;%raw
        mxste(d,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(d,i)   = c.myClasf.Zsc.fullSg.correct;%zscored
        mxzste(d,i)= c.myClasf.Zsc.fullSg.correctSTE;
        
        %get mean %correct for Nx-permuted dataset (chance)
        mxS(d,i) = c.myClasf.raw.fullSgShuf.correct; %raw
        mxzS(d,i) = c.myClasf.Zsc.fullSgShuf.correct; %zscored
        mxsteS(d,i) = c.myClasf.raw.fullSgShuf.correctSTE;
        mxzsteS(d,i) = c.myClasf.Zsc.fullSgShuf.correctSTE;
        
        %95% Confidence interval
        mxCI95S(d,i) = c.myClasf.raw.fullSgShuf.CI95updown(1);
        mxzCI95S(d,i) = c.myClasf.Zsc.fullSgShuf.CI95updown(1);
        
        %calculate p-value
        %statistical significance
        %of each % correct
        %pval is the percent of permutation accuracies
        %above data classification accuracy
        fprintf('%s \n','roi:',o.myROIname{i})
        fprintf('%s \n','---')
        
        %raw
        correctPerm = c.myClasf.raw.fullSgShuf.corrects;%raw
        pval.pval(d,i) = sum(correctPerm > mx(d,i))/length(correctPerm);
        fprintf('%s %.3f \n','p-value (raw):',pval.pval(d,i))
        
        %z-sc
        correctPermz = c.myClasf.Zsc.fullSgShuf.corrects;%z-sc
        pval.pvalz(d,i) = sum(correctPermz > mxz(d,i))/length(correctPermz);
        fprintf('%s %.3f \n','p-value (z-sc):',pval.pvalz(d,i))
        
        %notes
        pval.description = ['p-value for null hypothesis that classification',...
            '%correct comes from chance (N permutations)'];
        save('pval','pval')
    end
end

%----------------- %correct by coh and roi ---------------

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (red (95%CI))'])
    for co = 1 : nCoh
        hold all
        %chance
        bar(2*co-1,mxS(co,roi),'facecolor',[.5 .5 .5])
        %data
        bar(2*co,mx(co,roi),'facecolor',[0 0 0])
        %chance std and Confidence intervals
        %data classif acc set to 0
        %95%CI
        margErr = mxS(co,roi) - mxCI95S(co,roi);
        SLerrorbar([2*co-1 2*co],[mxS(co,roi) mx(co,roi)],...
            'yError',[margErr 0],'MarkerSize=1',...
            'Linestyle=none','Color=r')
        %STD
        %SLerrorbar([2*co-1 2*co],[mxS(co,roi) mx(co,roi)],...
        %    'yError',[mxsteS(co,roi) 0],'MarkerSize=1','Linestyle=none')
    end
    set(gca,'xtick',[1.5 3.5],'xticklabel',coi)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Coherence')
    hline(0.2);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')

%save figure
saveas(gcf, 'figureData', 'fig')

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (z-scored,red (95%CI))'])
    for co = 1 : nCoh
        hold all
        %chance
        bar(2*co-1,mxzS(co,roi),'facecolor',[.5 .5 .5])
        %data
        bar(2*co,mxz(co,roi),'facecolor',[0 0 0])
        %chance std and Confidence intervals
        %data classif acc set to 0
        %95%CI
        margErrz = mxzS(co,roi) - mxzCI95S(co,roi);
        SLerrorbar([2*co-1 2*co],[mxzS(co,roi) mxz(co,roi)],...
            'yError',[margErrz 0],'MarkerSize=1',...
            'Linestyle=none','Color=r')
        %std
        %SLerrorbar([2*co-1 2*co],[mxzS(co,roi) mxz(co,roi)],...
        %    'yError',[mxzsteS(co,roi) 0],'MarkerSize=1','Linestyle=none')
    end
    set(gca,'xtick',[1.5 3.5],'xticklabel',coi)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Coherence')
    hline(0.2);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')


%save figure
saveas(gcf, 'figureZsco', 'fig')

%table
for d = 1 : nCoh
    % %correct and ste by roi
    SLprintTable([mx(d,:)' mxste(d,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figAccurAndSTEbyRoi' num2str(d)], 'fig')
    % scored %correct and ste by roi
    SLprintTable([mxz(d,:)' mxzste(d,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figZscAccurAndSTEbyRoi' num2str(d)], 'fig')
    % pvals by roi
    SLprintTable([pval.pval(d,:)' pval.pvalz(d,:)'],rois,{'%c','ste'})
    saveas(gcf, ['pvals' num2str(d)], 'fig')
end
function [mx,mxz] = mySwitchaccAtTimeleaveOneOutStk(o)

rois = o.myROIname;

fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','Test p-value that %correct comes from chance (null Hyp)')
fprintf('%s \n','-------------------------------------------------------')

%class
myClass = [o.myClass{:}];

%path
analpath = [o.rootpath 'slStckAnalyses/' ...
    o.myGroup '/classif/' myClass '/' [o.myAnalysis{:}] '/'];

fprintf('%s \n','-------------')

%each roi
for i = 1 : o.nRois
    
    %check path
    if exist(analpath)
        %data path
        cd([analpath rois{i}]);
    else
        slPrintfStr('slfMRIclassifSAdisplay',' Path does not exist.')
    end
    
    %get file
    myfile = dir('ClassifStckSess*mat');
    clear('c')
    load(myfile.name)
    
    %store %correct and STE
    %actual data
    %Ndi by Nrois
    mx(i)    = c.myClasf.raw.fullSg.correct;%raw
    mxste(i) = c.myClasf.raw.fullSg.correctSTE;
    mxz(i)   = c.myClasf.Zsc.fullSg.correct;%zscored
    mxzste(i)= c.myClasf.Zsc.fullSg.correctSTE;
    
    %get mean %correct for Nx-permuted dataset (chance)
    mxS(i) = c.myClasf.raw.fullSgShuf.correct; %raw
    mxzS(i) = c.myClasf.Zsc.fullSgShuf.correct; %zscored
    mxsteS(i) = c.myClasf.raw.fullSgShuf.correctSTE;
    mxzsteS(i) = c.myClasf.Zsc.fullSgShuf.correctSTE;
    
    %95% Confidence interval
    mxCI95S(i) = c.myClasf.raw.fullSgShuf.CI95updown(1);
    mxzCI95S(i) = c.myClasf.Zsc.fullSgShuf.CI95updown(1);
    
    %calculate p-value
    %statistical significance
    %of each % correct
    %pval is the percent of permutation accuracies
    %above data classification accuracy
    fprintf('%s \n','roi:',o.myROIname{i})
    fprintf('%s \n','---')
    
    %raw
    correctPerm = c.myClasf.raw.fullSgShuf.corrects;%raw
    pval.pval(i) = sum(correctPerm > mx(i))/length(correctPerm);
    fprintf('%s %.3f \n','p-value (raw):',pval.pval(i))
    
    %z-sc
    correctPermz = c.myClasf.Zsc.fullSgShuf.corrects;%z-sc
    pval.pvalz(i) = sum(correctPermz > mxz(i))/length(correctPermz);
    fprintf('%s %.3f \n','p-value (z-sc):',pval.pvalz(i))
    
    %notes
    pval.description = ['p-value for null hypothesis that classification',...
        '%correct comes from chance (N permutations)'];
    save(['pval' o.myROIname{i}],'pval')
end

%----------------- %correct by coh and roi ---------------

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (red (95%CI))'])
    
    hold all
    %chance
    bar(1,mxS(roi),'facecolor',[.5 .5 .5])
    %data
    bar(2,mx(roi),'facecolor',[0 0 0])
    %chance std and Confidence intervals
    %data classif acc set to 0
    %95%CI
    margErr = mxS(roi) - mxCI95S(roi);
    SLerrorbar([1 2],[mxS(roi) mx(roi)],...
        'yError',[margErr 0],'MarkerSize=1',...
        'Linestyle=none','Color=r');
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    
end
legend('Chance(permutation)','data')
legend('boxoff')

%save figure
saveas(gcf, 'figAccData', 'fig')

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (z-scored,red (95%CI))'])
    
    hold all
    %chance
    bar(1,mxzS(roi),'facecolor',[.5 .5 .5])
    %data
    bar(2,mxz(roi),'facecolor',[0 0 0])
    %chance std and Confidence intervals
    %data classif acc set to 0
    %95%CI
    margErrz = mxzS(roi) - mxzCI95S(roi);
    SLerrorbar([1 2],[mxzS(roi) mxz(roi)],...
        'yError',[margErrz 0],'MarkerSize=1',...
        'Linestyle=none','Color=r');
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    hline(0.5);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')


%save figure
saveas(gcf, 'figAccZscored', 'fig')

%table
%correct and ste by roi
SLprintTable([mx(1,:)' mxste(1,:)'],rois,{'%c','ste'})
saveas(gcf, 'figAccurAndSTEbyRoi', 'eps')
% scored %correct and ste by roi
SLprintTable([mxz(1,:)' mxzste(1,:)'],rois,{'%c','ste'})
saveas(gcf,'figZscAccurAndSTEbyRoi', 'eps')
% pvals by roi
SLprintTable([pval.pval(1,:)' pval.pvalz(1,:)'],rois,{'%c','ste'})
saveas(gcf,'pvals','eps')
function [mx,mxz] = mySwitchbyDiraccAtTimeleaveOneOutStk(o)

rois = o.myROIname;
dire = [15 85 155 295];
ndire= length(dire);

%each coherence
fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','Test p-value that %correct comes from chance (null Hyp)')
fprintf('%s \n','-------------------------------------------------------')
for di = 1 : ndire
    %class
    myClass = ['mySwitch=1_x_myRandomDir=' num2str(dire(di)) 'mySwitch=2_x_myRandomDir=' num2str(dire(di))];
    %path
    analpath{di} = [o.rootpath 'slStckAnalyses/' ...
        o.myGroup '/classif/' myClass '/' [o.myAnalysis{:}] '/'];
    
    fprintf('%s %.2f \n','Dire:',dire(di))
    fprintf('%s \n','-------------')
    %each roi
    for i = 1 : o.nRois
        
        %check path
        if exist(analpath{di})
            %data path
            cd([analpath{di} rois{i}]);
        else
            slPrintfStr('slfMRIclassifSAdisplay',' Path does not exist.')
        end
        
        %get file
        myfile = dir('ClassifStckSess*mat');
        clear('c')
        load(myfile.name)
        
        %store %correct and STE
        %actual data
        %Ndi by Nrois
        mx(di,i)    = c.myClasf.raw.fullSg.correct;%raw
        mxste(di,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(di,i)   = c.myClasf.Zsc.fullSg.correct;%zscored
        mxzste(di,i)= c.myClasf.Zsc.fullSg.correctSTE;
        
        %calculate p-value
        %statistical significance
        %of each % correct
        %pval is the percent of permutation accuracies
        %above data classification accuracy
        fprintf('%s \n','roi:',o.myROIname{i})
        fprintf('%s \n','---')

    end
end

%----------------- %correct by dire and roi ---------------

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (red (95%CI))'])
    for di = 1 : ndire
        hold all
        %data
        bar(2*di,mx(di,roi),'facecolor',[0 0 0])        
    end
    set(gca,'xtick',[1.5 3.5 5.5 9.5],'xticklabel',dire)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Direction')
    hline(0.5);
    ylim([0 0.6])
end
legend('data')
legend('boxoff')

%save figure
saveas(gcf, 'figAccData', 'fig')

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (z-scored,red (95%CI))'])
    for di = 1 : ndire
        hold all
        %data
        bar(2*di,mxz(di,roi),'facecolor',[0 0 0])        
    end
    set(gca,'xtick',[1.5 3.5 5.5 9.5],'xticklabel',dire)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Direction')
    hline(0.5);
    ylim([0 0.6])
end
legend('data')
legend('boxoff')


%save figure
saveas(gcf, 'figAccZscored', 'fig')

%table
for di = 1 : ndire
    % %correct and ste by roi
    SLprintTable([mx(di,:)' mxste(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figAccurAndSTEbyRoi' num2str(di)], 'eps')
    % scored %correct and ste by roi
    SLprintTable([mxz(di,:)' mxzste(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figZscAccurAndSTEbyRoi' num2str(di)], 'eps')
    saveas(gcf, ['pvals' num2str(di)], 'eps')
end
function [mx,mxz] = mySwitchbyDirWithPermaccAtTimeleaveOneOutStk(o)

rois = o.myROIname;
dire = [15 85 155 295];
ndire= length(dire);

%each coherence
fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','Test p-value that %correct comes from chance (null Hyp)')
fprintf('%s \n','-------------------------------------------------------')
for di = 1 : ndire
    %class
    myClass = ['mySwitch=1_x_myRandomDir=' num2str(dire(di)) 'mySwitch=2_x_myRandomDir=' num2str(dire(di))];
    %path
    analpath{di} = [o.rootpath 'slStckAnalyses/' ...
        o.myGroup '/classif/' myClass '/' [o.myAnalysis{:}] '/'];
    
    fprintf('%s %.2f \n','Dire:',dire(di))
    fprintf('%s \n','-------------')
    %each roi
    for i = 1 : o.nRois
        
        %check path
        if exist(analpath{di})
            %data path
            cd([analpath{di} rois{i}]);
        else
            slPrintfStr('slfMRIclassifSAdisplay',' Path does not exist.')
        end
        
        %get file
        myfile = dir('ClassifStckSess*mat');
        clear('c')
        load(myfile.name)
        
        %store %correct and STE
        %actual data
        %Ndi by Nrois
        mx(di,i)    = c.myClasf.raw.fullSg.correct;%raw
        mxste(di,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(di,i)   = c.myClasf.Zsc.fullSg.correct;%zscored
        mxzste(di,i)= c.myClasf.Zsc.fullSg.correctSTE;
        
        %get mean %correct for Nx-permuted dataset (chance)
        mxS(di,i) = c.myClasf.raw.fullSgShuf.correct; %raw
        mxzS(di,i) = c.myClasf.Zsc.fullSgShuf.correct; %zscored
        mxsteS(di,i) = c.myClasf.raw.fullSgShuf.correctSTE;
        mxzsteS(di,i) = c.myClasf.Zsc.fullSgShuf.correctSTE;
        
        %95% Confidence interval
        mxCI95S(di,i) = c.myClasf.raw.fullSgShuf.CI95updown(1);
        mxzCI95S(di,i) = c.myClasf.Zsc.fullSgShuf.CI95updown(1);
        
        %calculate p-value
        %statistical significance
        %of each % correct
        %pval is the percent of permutation accuracies
        %above data classification accuracy
        fprintf('%s \n','roi:',o.myROIname{i})
        fprintf('%s \n','---')
        
        %raw
        correctPerm = c.myClasf.raw.fullSgShuf.corrects;%raw
        pval.pval(di,i) = sum(correctPerm > mx(di,i))/length(correctPerm);
        fprintf('%s %.3f \n','p-value (raw):',pval.pval(di,i))
        
        %z-sc
        correctPermz = c.myClasf.Zsc.fullSgShuf.corrects;%z-sc
        pval.pvalz(di,i) = sum(correctPermz > mxz(di,i))/length(correctPermz);
        fprintf('%s %.3f \n','p-value (z-sc):',pval.pvalz(di,i))
        
        %notes
        pval.description = ['p-value for null hypothesis that classification',...
            '%correct comes from chance (N permutations)'];
        save(['pval' o.myROIname{i}],'pval')
    end
end

%----------------- %correct by dire and roi ---------------

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (red (95%CI))'])
    for di = 1 : ndire
        hold all
        %chance
        bar(2*di-1,mxS(di,roi),'facecolor',[.5 .5 .5])
        %data
        bar(2*di,mx(di,roi),'facecolor',[0 0 0])
        %chance std and Confidence intervals
        %data classif acc set to 0
        %95%CI
        margErr = mxS(di,roi) - mxCI95S(di,roi);
        SLerrorbar([2*di-1 2*di],[mxS(di,roi) mx(di,roi)],...
            'yError',[margErr 0],'MarkerSize=1',...
            'Linestyle=none','Color=r');
    end
    set(gca,'xtick',[1.5 3.5 5.5 9.5],'xticklabel',dire)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Direction')
    hline(0.5);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')

%save figure
saveas(gcf, 'figAccData', 'fig')

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (z-scored,red (95%CI))'])
    for di = 1 : ndire
        hold all
        %chance
        bar(2*di-1,mxzS(di,roi),'facecolor',[.5 .5 .5])
        %data
        bar(2*di,mxz(di,roi),'facecolor',[0 0 0])
        %chance std and Confidence intervals
        %data classif acc set to 0
        %95%CI
        margErrz = mxzS(di,roi) - mxzCI95S(di,roi);
        SLerrorbar([2*di-1 2*di],[mxzS(di,roi) mxz(di,roi)],...
            'yError',[margErrz 0],'MarkerSize=1',...
            'Linestyle=none','Color=r');
        %std
        %SLerrorbar([2*co-1 2*co],[mxzS(co,roi) mxz(co,roi)],...
        %    'yError',[mxzsteS(co,roi) 0],'MarkerSize=1','Linestyle=none')
    end
    set(gca,'xtick',[1.5 3.5 5.5 9.5],'xticklabel',dire)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Direction')
    hline(0.5);
    ylim([0 0.6])
end
legend('Chance(permutation)','data')
legend('boxoff')


%save figure
saveas(gcf, 'figAccZscored', 'fig')

%table
for di = 1 : ndire
    % %correct and ste by roi
    SLprintTable([mx(di,:)' mxste(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figAccurAndSTEbyRoi' num2str(di)], 'eps')
    % scored %correct and ste by roi
    SLprintTable([mxz(di,:)' mxzste(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figZscAccurAndSTEbyRoi' num2str(di)], 'eps')
    % pvals by roi
    SLprintTable([pval.pval(di,:)' pval.pvalz(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['pvals' num2str(di)], 'eps')
end
function [mx,mxz] = mySwitchbyCohaccAtTimeleaveOneOutStk(o)

rois = o.myROIname;
dire = [0.06 0.12];
ndire= length(dire);

%each coherence
fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','Test p-value that %correct comes from chance (null Hyp)')
fprintf('%s \n','-------------------------------------------------------')
for di = 1 : ndire
    %class
    myClass = ['mySwitch=1_x_myRandomCoh=' num2str(dire(di)) 'mySwitch=2_x_myRandomCoh=' num2str(dire(di))];
    %path
    analpath{di} = [o.rootpath 'slStckAnalyses/' ...
        o.myGroup '/classif/' myClass '/' [o.myAnalysis{:}] '/'];
    
    fprintf('%s %.2f \n','Coh:',dire(di))
    fprintf('%s \n','-------------')
    %each roi
    for i = 1 : o.nRois
        
        %check path
        if exist(analpath{di})
            %data path
            cd([analpath{di} rois{i}]);
        else
            slPrintfStr('slfMRIclassifSAdisplay',' Path does not exist.')
        end
        
        %get file
        myfile = dir('ClassifStckSess*mat');
        clear('c')
        load(myfile.name)
        
        %store %correct and STE
        %actual data
        %Ndi by Nrois
        mx(di,i)    = c.myClasf.raw.fullSg.correct;%raw
        mxste(di,i) = c.myClasf.raw.fullSg.correctSTE;
        mxz(di,i)   = c.myClasf.Zsc.fullSg.correct;%zscored
        mxzste(di,i)= c.myClasf.Zsc.fullSg.correctSTE;
        
        fprintf('%s \n','roi:',o.myROIname{i})
        fprintf('%s \n','---')

    end
end

%----------------- %correct by dire and roi ---------------

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi}])
    for di = 1 : ndire
        hold all                
        SLerrorbar(2*di,mx(di,roi),...
            'yError',mxste,'MarkerSize=1',...
            'Linestyle=none','Color=r');
                %data
        bar(2*di,mx(di,roi),'facecolor',[0 0 0])
    end
    set(gca,'xtick',[1.5 3.5],'xticklabel',dire)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Coherence')
    hline(0.5);
    ylim([0 1])
end
legend('data')
legend('boxoff')

%save figure
saveas(gcf, 'figAccData', 'fig')

figure('color','w')
for roi = 1 : o.nRois
    subplot(1,o.nRois,roi)
    title([o.myROIname{roi} ' (z-scored)'])
    for di = 1 : ndire        
        hold all                
        SLerrorbar(2*di,mxz(di,roi),...
            'yError',mxzste,'MarkerSize=1',...
            'Linestyle=none','Color=r');
                %data
        bar(2*di,mxz(di,roi),'facecolor',[0 0 0])
        %std
        %SLerrorbar([2*co-1 2*co],[mxzS(co,roi) mxz(co,roi)],...
        %    'yError',[mxzsteS(co,roi) 0],'MarkerSize=1','Linestyle=none')
    end
    set(gca,'xtick',[1.5 3.5],'xticklabel',dire)
    
    %label
    if roi == 1
        ylabel('Accuracy (%correct)')
    end
    xlabel('Coherence')
    hline(0.5);
    ylim([0 1])
end
legend('data')
legend('boxoff')


%save figure
saveas(gcf, 'figAccZscored', 'fig')

%table
for di = 1 : ndire
    % %correct and ste by roi
    SLprintTable([mx(di,:)' mxste(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figAccurAndSTEbyRoi' num2str(di)], 'eps')
    % scored %correct and ste by roi
    SLprintTable([mxz(di,:)' mxzste(di,:)'],rois,{'%c','ste'})
    saveas(gcf, ['figZscAccurAndSTEbyRoi' num2str(di)], 'eps')    
end