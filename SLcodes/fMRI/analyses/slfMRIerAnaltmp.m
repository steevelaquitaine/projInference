
%slfMRIerAnaltmp.m
%
%
%     Author: Steeve Laquitaine
%       date: 150909
%    purpose: run some basic event-related analyses
%
%      usage:
%
%               o.myVar = 'myRandomCoh';
%           o.myROIname = {'rV3A'};
%             o.myGroup = 'MotionComp';
%      o.myGroupScanNum = 1;
%            o.scanList = 1:18;         
%          o.myr2Thresh = 0.2;
%             o.taskNum = 2;
%            o.phaseNum = 1;       
%          o.segmentNum = 2;                    
%              o.hdrlen = 50;                   
%              o.myBase = 's0025_flatR_WM_occipital_Rad90.hdr';
%             o.myPath = '~/dataold/datafMRI/sltaskdotdirfmri05/s02520150814';
%          o.myBasePath = '~/dataold/datafMRI/mlrAnatDB/s0025/mlrBaseAnatomies/';
%         o.myCanonical = '~/dataold/datafMRI/mlrAnatDB/s0025/surfaces/s0025_mprage_pp.nii';
%           o.myROIpath = '~/dataold/datafMRI/mlrAnatDB/s0025/mlrROIs/';
%         
%         o = slfMRIerAnaltmp(o)
%
%
%  Analyses:
%
%   - map r2 calculated from fitting deconvolved hdrs with Bold time Series
%     (% signal changes on a flat map) for all or separated conditions
%   - bar # voxels by r2 threshold.
%   - plot er by voxels sorted by r2.
%
%
%  Reference : 
%
%   http://gru.stanford.edu/doku.php/mrTools/scriptingExamplesEventRelated

function [o,v] = slfMRIerAnaltmp(o)

tic

%backup m.file with date and time
%filename = matlab.desktop.editor.getActiveFilename;
%SLBackup(filename)

%launch mrLoadRet
cd(o.myPath)
mrLoadRet([]);
v = getMLRView;

%------------------------------------
%        Load Base anatomy 
%------------------------------------
BaseFileType = o.myBase(find(o.myBase=='.') : end);

%case flat
if strcmp(BaseFileType,'.hdr')
    v = loadAnat(v,o.myBase,o.myBasePath);
end

%case surface
if strcmp(BaseFileType,'.off')
    base = slImportSurfaceOFF([o.myBasePath o.myBase],[],o.myCanonical);
    v = viewSet(v,'newbase', base);
end
groupNames = viewGet(v,'groupNames');            %check that concatenation exists
refreshMLRDisplay(v);                            %refresh view


%--------------------------------  er ---------------------------
v = viewSet(v,'curGroup',o.myGroup);
v = viewSet(v,'curScan',o.myGroupScanNum);

%check if er. exists
fprintf('\n %s \n \n','-------------------- README -------------------')
YorN = input(['!! README !! er. takes a while for long scans (e.g., concatenation)',...
    '. If you"ve run er LOAD IT !! Otherwise make sure that you set',...
    'maxBlocksize in Edit - preferences to a reasonable low value',...
    '(e.g., 30GB) Do you want to load an already existing "erAnal" ',...
    'for speed ? [y/n] :'],'s');

%load if exists
if strcmp(YorN,'y')
    fprintf('\n %s \n \n','------------ Loading existing erAnal ----------')
    v = loadAnalysis(v);
else
    %or run er
    fprintf('\n %s \n \n',' ------------ Running er analysis ------------ ')
    
    %init deflt params
    [~,params] = eventRelated(v,[],'justGetParams=1','defaultParams=1',...
        'scanList',o.myGroupScanNum);                         
    fprintf('%s %s %s \n','(slfMRIerAnal)','Available variables for er analysis are :',params.scanParams{1}.paramInfo{7}{3})
    params.applyFiltering                                          = 1;           %filtering
    params.scanParams{o.myGroupScanNum}.stimvolVarInfo.varname     = o.myVar;     %er variable
    params.scanParams{o.myGroupScanNum}.stimvolVarInfo.taskNum     = o.taskNum;   %task
    params.scanParams{o.myGroupScanNum}.stimvolVarInfo.segmentNum  = o.segmentNum;
    %params.scanParams{o.myGroupScanNum}.varname                    = o.myVar;     %er variable
    %params.scanParams{o.myGroupScanNum}.taskNum                    = o.taskNum;   %task
    %params.scanParams{o.myGroupScanNum}.segmentNum                 = o.segmentNum;
    [v,d] = eventRelated(v,params); %r2 map
end

%-------------------------- Display all ROIs -----------------------------
fprintf('\n %s \n \n','----------------- Loading all Rois ----------------')
slLoadAllROI(v,o.myROIpath,1);


%-----------------------  PARFOR plots over ROIS --------------------------
o.erCav      = [0 0 0; .7 .7 .7];
o.scanDims   = viewGet(v,'scanDims'); %N by N by # slice
o.r2All      = viewGet(v,'overlayData',o.myGroupScanNum);%r2 for all vox in volume
o.nRois      = length(o.myROIname);

%stimulus convolution matrix (Nvols by NhdrFrames with 1 and 0)
[stimvol,o.stimNames] = getStimvol(v,o.myVar,'taskNum',o.taskNum,...
    'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);                     
o.nhdr = length(stimvol);
o.scm  = makescm(v,o.hdrlen,1,stimvol);

%get task segments for plot
a = viewGet(v,'stimfile');

%mean seg. duration
for i = 1: length(a)
    o.sglen(i,:) = mean([a{i}.task{o.taskNum}{1}.segmin; a{i}.task{o.taskNum}{1}.segmax]); 
end
o.sglen = unique(o.sglen,'rows');

%check that unique duration per segment was found
if size(o.sglen,1) ~= 1
    fprintf('%s \n','(slfMRIerAnaltmp) The scans trials have different',...
        'segment lengths. e.g., trial 1 : [1 0.3 5 3 11.7] and trial N : ',...
        '[1 3 10 3 20]. The code currently requires all scans to have ',...
        'same segment lengths across trials. e.g., all trials :',...
        '[1 0.3 5 3 11.7]')
    keyboard
end
cssglen  = [0 cumsum(o.sglen)];
o.TrlEvT = cssglen - cssglen(o.segmentNum);
evcol    = linspecer(length(o.sglen));

%------------------------- load ROITseries -----------------------------------
%This is insanely slow, lots of memory !! make sure edit - preferences - maxBlocksize is low enough
%GET rid of it asap. Prob is er analysis doesn't load and store the
%entire volume but load slice by slice. Thus cannot be used here to get
%the ROI TSeries. The fastest way so far is to save the ROI TSeries once
%run once and reload them instead of loadROITSeries them.
%e.myROI.scanCoords = getROICoordinates(v,myROIname);
YesorNo = input(['!! README !! Loading ROI tseries takes a while for long'...
    ' scans (e.g., concatenation). Do you want to load existing ROItseries ? [y/n]'],'s');
%if no load roi tseries
if strcmp(YesorNo,'n')
    for ROIi = 1 : o.nRois
        myROIi = [];
        fullroi = [o.myROIpath o.myROIname{ROIi} '.mat']; %roi with path
        myROIi = loadROITSeries(v,fullroi);        %load (very slow)
        mkdir(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])    %save once to reload
        cd(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])
        save(['tSeries' o.myGroup o.myROIname{ROIi} date],'myROIi')
        cd ..
        e.myROIi{ROIi} = myROIi; %prep for analyses
    end
end
%otherwise load roi tseries
if strcmp(YesorNo,'y')
    for ROIi = 1 : o.nRois
        e.myROIi{ROIi} = importdata(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/tSeries' o.myGroup o.myROIname{ROIi} date '.mat']);
    end
end

for ROIi = 1 : o.nRois
    %cleanup. Keep view, data, inputs
    keep = {'v','e','o','ROIi','keep','a','i'}; a = whos;
    for i = 1 : length(a)
        if ~sum(strcmp(a(i).name,keep))>0
            clear(a(i).name)
        end
    end
    clear('keep','a','i')
    myROIname  = o.myROIname{ROIi};
    figure;
    f = gcf;
    set(gcf,'color','w')
    
    %plot (er average over r2>thresh voxels)
    %---------------------------------------
    %make an er plot from ROI r2's thresholded map. We keep the
    %time series of the voxels with high r2 (vox. activity
    %that are modeled best). We average those timeSeries into on time series
    %and we use this time series to estimate an hemodynamic response function
    %for our different stimli (e.g., coherences).
    %make an average time Series for all the voxels that have a r2 greater than
    %a certain threshold.
    subplot(1,4,1)
    ercolors = [.5 0 0; 1 0 0];
    e.myROIi{ROIi}.linCoor         = sub2ind(o.scanDims,e.myROIi{ROIi}.scanCoords(1,:),e.myROIi{ROIi}.scanCoords(2,:),e.myROIi{ROIi}.scanCoords(3,:)); %ROI coord converted to linear
    e.myROIi{ROIi}.r2All           = o.r2All(e.myROIi{ROIi}.linCoor); %pick ROI's r2s
    e.tS                    = mean(e.myROIi{ROIi}.tSeries(e.myROIi{ROIi}.r2All>o.myr2Thresh,:),1); %mean TSeries for high r2 vox (Nvols by 1)
    e.nVoxAll               = sum(e.myROIi{ROIi}.r2All>0); %# voxels with r2 greater than the threshold
    e.nVoxAbThres           = sum(e.myROIi{ROIi}.r2All > o.myr2Thresh);
    e.d                     = getr2timecourse(e.tS,o.nhdr,o.hdrlen,o.scm,viewGet(v,'framePeriod'));                                          %deconv. hdrs (fitTimecourse)
    for j = 1 : o.nhdr
        SLerrorbar(e.d.time, e.d.ehdr(j,:),'yError', e.d.ehdrste(j,:),...
            ['Color=[',num2str(ercolors(j,:)),']'],'MarkerSize=7','linewidth',2,...
            'linesmoothing','on');
    end
    hline(0)    
    xlim([min(o.TrlEvT)-1 max([o.TrlEvT e.d.time(end)])])
    xlabel('Time (sec)')
    ylabel('Bold (% signal change)')
    
    %add tasks segments
    for i = 1 : length(o.sglen)+1
        vline(o.TrlEvT(i))
    end
    h = get(gca,'children');
    if ~strcmp(o.stimNames,'_all_')
        legend(h((o.hdrlen+1)*[1:o.nhdr]),o.stimNames)
    else
        legend('High r2 voxels')
    end
    
    
    %plot (er averg over all voxels)
    %------------------------------
    subplot(1,4,2)
    e.tS           = mean(e.myROIi{ROIi}.tSeries);                                    %average TSeries
    e.dall         = getr2timecourse(e.tS,o.nhdr,o.hdrlen,o.scm,viewGet(v,'framePeriod')); %deconv. hdrs (fitTimecourse)
    for j = 1 : o.nhdr
        SLerrorbar(e.dall.time, e.dall.ehdr(j,:), 'yError', e.dall.ehdrste(j,:),...
            ['Color=[',num2str(o.erCav(j,:)),']'],'MarkerSize=7','linewidth',2,...
            'linesmoothing','on');
    end
    
    hline(0)    
    box off
    set(gca,'fontsize',14)
    title(sprintf('%s %s %s %i %s \n %s %.2g %s %i %s','er:',myROIname,', hdrlen=',o.hdrlen,...
        'vols','all voxels(gray)'),'fontsize',10)
    ylabel('Bold (% signal change)','fontsize',10)
    xlabel('Time (seconds)','fontsize',10)
    xlim([min(o.TrlEvT)-1 max([o.TrlEvT e.d.time(end)])])
    set(gca,'fontsize',10)
    h = get(gca,'children');
    if ~strcmp(o.stimNames,'_all_')
        legend(h((o.hdrlen+1)*[1:o.nhdr]),o.stimNames)
    else
        legend(h([o.hdrlen+1 2*(o.hdrlen+1)]),'All voxels')
    end
    ylim([-4 4])
    vline(0)
    
    
    %plot r2 density
    %---------------
    subplot(1,4,3)
    [n,x] = hist(e.myROIi{ROIi}.r2All,0:0.1:1);
    SLdrawBar(n/e.nVoxAll *100,x,1:11,'text',n,'FaceColor',[.5 .5 .5])
    title(sprintf('%s %.2f %s %.2f %s \n %i %s',' Mean r2 : ',mean(e.myROIi{ROIi}.r2All),...
        '(',std(e.myROIi{ROIi}.r2All),'std)',e.nVoxAll ,' vox.'))
    xlabel('r2','fontsize',10)
    ylabel('Occurence (% of voxels)','fontsize',10)
    set(gca,'fontsize',10)
    ylim([0 100])
    
    
    
    %plot # of voxels by r2 threshold
    %--------------------------------
    subplot(1,4,4)
    r2thres = 0:0.1:0.8;
    e.nVoxByR2thres = nan(length(r2thres),1); %# of voxels
    for j = 1 : length(r2thres)
        e.nVoxByR2thres(j) = sum(e.myROIi{ROIi}.r2All>r2thres(j));
    end
    SLdrawBar(e.nVoxByR2thres, r2thres, 1 : length(r2thres))
    xlabel('r2 thresholds','fontsize',10)
    ylabel('# of voxels with r2 > threshold','fontsize',10)
    title(sprintf('%i %s %s \n %.0f %s %.2f %s',e.nVoxAll ,' voxels',myROIname,e.nVoxAbThres/e.nVoxAll *100,'% vox (r2>',o.myr2Thresh,')'),'fontsize',10)
    set(gca,'fontsize',10)
    
    
    %----------------
    %plot er by voxel
    %----------------
    e.ix                     = 1 : numel(e.myROIi{ROIi}.r2All);%sort voxels by r2
    [~,e.ix]                 = sortrows(e.myROIi{ROIi}.r2All');
    e.ix                     = e.ix(end:-1:1);
    e.myROIi{ROIi}.tSr2Sorted       = e.myROIi{ROIi}.tSeries(e.ix,:);   %sort TSeries
    C                        = linspecer(e.nVoxAll );
    
%     %-------------------
%     %Plot each voxel hdr
%     %-------------------
%     ercolors = [.5 0 0; 1 0 0];
%     nVoxperfig = 100;
%     myAx = repmat(1:nVoxperfig,1,ceil(e.nVoxAll /nVoxperfig));
%     
%     %each voxel
%     for Voxi = 1 : e.nVoxAll 
%         
%         %new fig. every 100 vox
%         if myAx(Voxi) == 1
%             figure('color','w')
%         end
%         subplot(sqrt(nVoxperfig),sqrt(nVoxperfig),myAx(Voxi));
%         tSVoxi= e.myROIi{ROIi}.tSr2Sorted(Voxi,:);
%         dVoxi = getr2timecourse(tSVoxi,o.nhdr,hdrlen,o.scm,viewGet(v,'framePeriod'),0);
%         for j = 1 : o.nhdr
%             SLerrorbar(dVoxi.time,dVoxi.ehdr(j,:),'yError',dVoxi.ehdrste(j,:),...
%                 ['Color=[',num2str(ercolors(j,:)),']'],...
%                 'MarkerSize=5',...
%                 'linewidth',2,...
%                 'linesmoothing','on');
%             xlim([0 dVoxi(1).time(end)])
%             ylim([-3 +3])
%         end
%         hold on; plot([0 dVoxi(1).time(end)],[0 0],'color',[.5 .5 .5])
%         annotation('textbox', [0 0.9 1 0.1], ...
%             'String', [myROIname{ROIi},'(', e.nVoxAll ,')'], ...
%             'EdgeColor', 'none', ...
%             'HorizontalAlignment', 'center','fontsize',10)
%         
%         title(['r2:',num2str(e.myROIi{ROIi}.r2All(e.ix(Voxi)))])
%         xlabel('Time (sec)','fontsize',10)
%         ylabel('% Signal change','fontsize',10)
%         drawnow
%     end
    
    o.duration = toc;
    save(['erBy' o.myVar date])
    
    %backup
    fprintf('\n')
    fprintf('%s \n','------------------Saving data---------------')
    mkdir(['myErAnal' o.myROIname{ROIi}])
    cd(['myErAnal' o.myROIname{ROIi}])
    save(['myErAnal' o.myVar o.myROIname{ROIi} date],'o','e')
    saveas(f, ['myErAnal' o.myVar o.myROIname{ROIi}], 'fig')
    o.duration = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP over ROIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









