
%slfMRIerAnal.m
%
%
%     Author: Steeve Laquitaine
%       date: 150909
%    purpose: run some basic event-related analyses
%
%      usage: 
%
%         %params  
%         o.myVar         = '_all_'; 
%         o.taskNum       = 2;
%         o.phaseNum      = 1;
%         o.segmentNum    = 2;
%         o.myPath        = '~/data/sltaskdotdirfmri05/s02520150814';
%         o.myBase        = 's0025_flatR_WM_occipital_Rad90.hdr';
%         o.myBasePath    = '~/data/mlrAnatDB/s0025/mlrBaseAnatomies/';
%         o.myROIname     = {'rMT'}%,'lhV4','lIPS1','lIPS2','lIPS3','lV1',...
%                           %'lV2d','lV2v','lV3A','lV3d','lV3v','lV7',...
%                           %'rhV4','rIPS1','rIPS2','rIPS3','rLO1','rLO2',...
%                           %'rMT','rV1','rV2d','rV2v','rV3A','rV3B','rV3d',...
%                           %'rV3v','rV7'};
%         o.myROIpath     = '~/data/mlrAnatDB/s0025/mlrROIs/'
%         o.myGroup       = 'MotionComp'; %'Concatenation';
%         o.scanList      = 1 : 18;
%         o.myGroupScanNum = 1;
%         o.myr2Thresh     = 0.2;
%         o.hdrlen       = 50; %(in frames)

%         %run
%         o = slfMRIerAnal(o)
% 
%
%  Analyses:
%
%   - map r2 calculated from fitting deconvolved hdrs with Bold time Series 
%     (% signal changes on a flat map)
%   
%   - bar # voxels by r2 threshold.
%
%   - plot er by voxels sorted by r2.
%
%
%--------
%Analyses
%--------
%
%
%
%------
%To do: 
%------
%
% - Xlabel, ylabels,title,legend,summary info of the preprocessing and analysis
%   done, explanation errorbar.
% - add task event markers on plots (motion onset, response, ITI)
% - r2 map on full brain surface
% - r2 map on occipt + parietal flat 
% - total r2 density and r2 basic stats, mean, std, skew (all brain)
% - repeat per ROI
% - Occiptal average hdr - average per roi - each vox.
% - r2 statistics (density, mean,std)
% - Parietal average hdr - same
% - r2 statistics (density, mean,std)

%- load surface and flat map for each hemisphere
%- put two figures (one for each hemisphere) with subplots (analyses) with
%other figures of individual voxels edr for each hemisphere and roi

 
function o = slfMRIerAnal(o)

tic 

%backup m.file with date and time
%filename = matlab.desktop.editor.getActiveFilename;
%SLBackup(filename)

%launch mrLoadRet
cd(o.myPath)
mrLoadRet([])                           
v = getMLRView;

%------------------------------------
%%%%%% Load Base anatomy %%%%%%%%%%%%
%------------------------------------
BaseFileType= o.myBase(find(o.myBase=='.') : end);

%case flat
%---------
if strcmp(BaseFileType,'.hdr')
    v = loadAnat(v,o.myBase,o.myBasePath); 
end

%case surface
%------------
if strcmp(BaseFileType,'.off')
    
    %base = gruImportSurfaceOFF([o.myBasePath o.myBase],'bothHemiFlag',o.myCanonical); 
    base = slImportSurfaceOFF([o.myBasePath o.myBase],[],o.myCanonical); 
    v = viewSet(v,'newbase', base);
    
end
groupNames = viewGet(v,'groupNames');            %check that concatenation exists
refreshMLRDisplay(v);                            %refresh view

%---------------------------------------
%%%%%% run er analysis %%%%%%%%%%%%%%%%%%   
%---------------------------------------

%make sure we have concat
if sum(strcmp(groupNames,'Concatenation'))~=1;

    %set motion compensation group for concatenation
    curGroup = 'MotionComp';
    v = viewSet(v,'curGroup',curGroup);
    
    %- detrend, high-pass filter 
    %- get rid of junk frames, 
    %- do a de-trending step using a hi-pass filter if you
    %specify that, and then concatenate the resulting time series into a
    %single time series that will be saved into a new group called
    %Concatenation. The mat files are saved in folder "concatenation".
    %set concatenation parameters
    [~,params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1','scanList',o.scanList);
    
    %concatenate (stack preprocessed scans into one long scan)
    concatTSeries(v,params);
    fprintf('\n %12s \n','(slfMRIerAnal) Motion compensation scans have been concatenated')    
end
v = viewSet(v,'curGroup',o.myGroup);
v = viewSet(v,'curScan',o.myGroupScanNum);

%er
%--
[~,params] = eventRelated(v,[],'justGetParams=1','defaultParams=1','scanList',o.myGroupScanNum);                         %init default params
fprintf('%s %s %s \n','(slfMRIerAnal)','Available variables for er analysis are :',params.scanParams{1}.paramInfo{7}{3})
params.scanParams{o.myGroupScanNum}.stimvolVarInfo.varname     = o.myVar;  %er variable
params.applyFiltering                                          = 1;        %filtering
params.scanParams{o.myGroupScanNum}.stimvolVarInfo.taskNum     = o.taskNum;%task
params.scanParams{o.myGroupScanNum}.stimvolVarInfo.segmentNum  = o.segmentNum;
v = eventRelated(v,params);                                       %r2 map

%------------
%load all ROI
%------------
v   = slLoadAllROI(v,o.myROIpath,1);
img = refreshMLRDisplay(viewGet(v,'viewNum')); %separate figure
figure;
image(img);
axis square;
axis off;

%-----
%PLOTS
%-----
%---------------------------------------
%plot (er average over r2>thresh voxels)
%---------------------------------------
%make an er plot from ROI r2's thresholded map. We keep the
%time series of the voxels with high r2 (vox. activity 
%that are modeled best). We average those timeSeries into on time series
%and we use this time series to estimate an hemodynamic response function
%for our different stimuli (e.g., coherences).
%make an average time Series for all the voxels that have a r2 greater than
%a certain threshold.
myROI          = loadROITSeries(v,o.myROIname);
scanDims       = viewGet(v,'scanDims');             %N by N by # slice 
myROI.linCoor  = sub2ind(scanDims,myROI.scanCoords(1,:),myROI.scanCoords(2,:),myROI.scanCoords(3,:)); %ROI coord converted to linear
r2             = viewGet(v,'overlayData',o.myGroupScanNum);    %r2 value for every vox in the volume (Nheight by Nwidth by Nslices)
myROI.r2       = r2(myROI.linCoor);                 %pick ROI's r2s
tS             = mean(myROI.tSeries(myROI.r2>o.myr2Thresh,:)); %average TSeries for high r2 vox (Nvols by 1)
o.nVoxAll      = sum(myROI.r2>0);                   %# voxels with r2 greater than the threshold
o.nVoxAbThres  = sum(myROI.r2 > o.myr2Thresh);
stimvol        = getStimvol(v,o.myVar,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum); %Ntrials by 1 hdr from the average TSeries.
nhdr           = length(stimvol);
scm            = makescm(v,o.hdrlen,1,stimvol);     %make stim conv. matrnVoxAbThres (Nvols by NhdrFrames with 1 and 0)
o.d            = getr2timecourse(tS,nhdr,o.hdrlen,scm,viewGet(v,'framePeriod')); %deconv. hdrs (fitTimecourse)

figure('color','w');
hold all
ercolors = [.9 0 0; 1 0 0];
for j = 1 : nhdr
    SLerrorbar(o.d.time, o.d.ehdr(j,:), 'yError', o.d.ehdrste(j,:),['Color=[',num2str(ercolors(j,:)),']'],'linewidth',2,'linesmoothing','on')
end
hold on; plot([0 o.d.time(end)],[0 0],'color',[.5 .5 .5])
box off
set(gca,'fontsize',14)
title(sprintf('%s %i %s %s %s %2g %s %i %s','Event-related:',o.nVoxAbThres,'voxels ',o.myROIname{:},'r2>',o.myr2Thresh,', hdrlen=',o.hdrlen,'vols'))
ylabel('% Signal change')
xlabel('Time (seconds)')

%------------------------------------
%plot # of voxels with r2 > r2 thresh
%------------------------------------
r2thres = 0:0.05:0.8;
nVoxByR2thres = nan(length(r2thres),1); %# of voxels
for j = 1 : length(r2thres)
    nVoxByR2thres(j) = sum(myROI.r2>r2thres(j));
end

%draw the # of voxels against threshold
figure
SLdrawBar(nVoxByR2thres, r2thres, 1 : length(r2thres))
xlabel('r2 thresholds')
ylabel('# of voxels with r2 > threshold')
title([num2str(o.nVoxAll),'voxels, ',o.myROIname])

%----------------
%plot er by voxel
%----------------
ix = 1 : numel(myROI.r2);              %sort voxels by r2
[~,ix] = sortrows(myROI.r2');
ix = ix(end:-1:1);
myROI.tSr2Sorted = myROI.tSeries(ix,:);%sort TSeries

%figure
myscreensz = get(0,'screensize');
figure('color','w','position',[0.25 0.25 0.25 1].*myscreensz([3 4 3 4]));

ehdrVals = [];

%each voxel
for Voxi = 1 : o.nVoxAll
    
    h(Voxi)       = subplot(ceil(sqrt(o.nVoxAll)),ceil(sqrt(o.nVoxAll)),Voxi);
    tSVoxi        = myROI.tSr2Sorted(Voxi,:);
    o.dVoxi(Voxi) = getr2timecourse(tSVoxi,nhdr,o.hdrlen,scm,viewGet(v,'framePeriod')); %deconvolve the voxel's hdr

    %each stim condition (e.g., coherence)
    for j = 1 : nhdr

        SLerrorbar(o.dVoxi(Voxi).time,o.dVoxi(Voxi).ehdr(j,:),'yError',o.dVoxi(Voxi).ehdrste(j,:),...
            ['Color=[',num2str(ercolors(j,:)),']'],...
            'linewidth',2,...
            'linesmoothing','on',...
            'MarkerSize',5)

        xlim([0 max(o.dVoxi(Voxi).time)])
        ylim([min(o.dVoxi(Voxi).ehdr(:)) max(o.dVoxi(Voxi).ehdr(:))])

    end 
    
    hold on; plot([0 o.dVoxi(1).time(end)],[0 0],'color',[.5 .5 .5])
    
    %axis square
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', [o.myROIname{:},'(', o.nVoxAll,')'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

    title([' - r2:',num2str(myROI.r2(ix(Voxi)))])
    if Voxi==o.nVoxAll
        xlabel('Time (sec)')
        
    end
    if Voxi==1
        ylabel('% Signal change')
    end
    %max for plots
    ehdrVals = [ehdrVals; o.dVoxi(Voxi).ehdr(:)];
    
end
o.duration = toc;
set(h,'xlim',[0 o.dVoxi(1).time(end)])
set(h,'ylim',[min(ehdrVals) max(ehdrVals)])

v = refreshMLRDisplay(v);

%arrange fig on screen
SLpositionFigures










