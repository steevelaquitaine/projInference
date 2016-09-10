
%Author: Steeve Laquitaine
%date: 140428 last modification 140501 13:08
%file name: SLfMRIcorrelationAnalysis

%usage: The correlation analysis produces a coherence map on a flat
%map. The analysis is fully scripted. You can just run it.
%It does not requires the GUI except at one point.

%purpose: Run correlation analysis on the localizer experiment. This
%analysis identifies motion responsive areas as the areas showing the
%highest correlation with the alternation between stimulus - no
%stimulus period of the task (in our case: "coherent-random motion" and
                              %"coherent motion - black screen").

%%backup command window
%diary Log_SLfMRIcorrelationAnalysis_140501.txt

%%backup m.file with date and time
%filename = matlab.desktop.editor.getActiveFilename;
%SLBackup(filename)

%Run correlaton analysis on the localizer experiment.

%path
%cd /Users/steeve/fMRI_data/s02520140501/
cd ~/data/sltaskdotdirfmri05/s02520150814

%launch view
%mrLoadRet([])
v = getMLRView;

%set group to analyse (notched filtered motion compensated data)
%curGroup = 'MotionComp';
curGroup               = 'MotionComp';
%viewSet(v,'curGroup',curGroup);
%scanNumFromDescription = 'randMotion';    %scan to analyse (e.g., localizer with random motion, scaNum 1)
scanNumFromDescription = 'McompFullMotionLoc01';  %scan to analyse (e.g., localizer with random motion, scaNum 1)
%scanNumFromDescription='black';
%scanNum = viewGet(v,'scanNumFromDescription',scanNumFromDescription,curGroup);
%viewSet(v,'curScan',scanNum);
o.myFlat               = 's0025_flatR_WM_occipital_Rad90.hdr'; %base anatomy
o.myFlatPath           = '~/data/mlrAnatDB/s0025/mlrBaseAnatomies/';

%set views
scanNum = viewGet(v,'scanNumFromDescription',scanNumFromDescription,curGroup);
viewSet(v,'curGroup',curGroup);
viewSet(v,'curScan',scanNum);

%Init analysis
[v,params] = corAnal(v,[],'justGetParams=1','defaultParams=1');

%params analysis
params.recompute(scanNum) = 1;
params.ncycles(scanNum)   = 11;

%----
%run
%----
corAnal(v,params);

keyboard

%load stuffs
loadAnalysis(v,'corAnal/corAnal.mat'); %load analysis
loadAnat(v,o.myFlat,o.myFlatPath);     %load flat map anatomy 

%get and set coherence overlay (coh val by voxel)
overlayNames = viewGet(v,'overlayNames');       
OverlayNum   = find(strcmp('co',overlayNames));
viewSet(v,'curOverlay',OverlayNum);
analysisNum  = viewGet(v,'currentAnalysis'); %current analysis (e.g., corAnal)

%get data to overlay (e.g., coherence)
overlaydata = viewGet(v,'overlaydata',scanNum,OverlayNum,analysisNum);

%Display overlay (map appears)
%get group
groupNum = viewGet(v,'currentGroup');

%display
MyOverlayName = ['overlayName=My',cell2mat(overlayNames(OverlayNum))];
mrDispOverlay(overlaydata,scanNum,groupNum,v,MyOverlayName,...
              'colormapType','normal',...
              'cmap',hot(256));

%Load and overlay ROIs
%right
loadROI(v,{'rhV4','rV3v','rV2v','rV1','rV2d','rV3d',...
    'rV3a','rV3b','rLO1','rLO2','rMT'},0,...
        '/Users/steeve/fMRI_data/s02520140220/ROIs')
%left
%loadROI(v,{'lhV4','lV3v','lV2v','lV1','lV2d','lV3d',...
    %     'lV3a','lV3b','lLO1','lLO2','lMT'},0,...
          %     '/Users/steeve/fMRI_data/s02520140220/ROIs')

%show ROIs
viewSet(v,'showrois','all perimeter');

%color ROIs
roiNum=viewGet(v,'roiGroup');
for j = 1 : numel(roiNum)
    viewSet(v,'roiColor','white',roiNum(j));
end

%display everything altogether
img=refreshMLRDisplay(viewGet(v,'viewNum'));

%draw results in a matlab figure
% figure;
% image(img);
% axis square;
% axis off;

%backup and end diary
%save([filename(1:end-2),'.mat'])
%diary off




