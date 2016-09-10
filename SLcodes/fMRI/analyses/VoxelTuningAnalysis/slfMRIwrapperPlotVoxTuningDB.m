
%slfMRIwrapperPlotVoxTuningDB.m
%
%
% author: steeve laquitaine
%   date: 160310
%purpose: plot voxel tuning visual database (e.g., to stimulus motion direction)
%         e.g., stimulus motion direction tuning curves for each voxel
%
%Description: 
%   1) initialize data and analysis info
%   2) make database
%   3) vizualize sorted voxel tunings by switch variable

%% initialize analysis and make fMRI database
slfmriInitAnalysisTaskDotDirfMRI05
[d,o] = slfMRImakeDatabase(o,varargin);

%% checks
%%% simulate data with equiprobable distribution of direction preferences
%d = slsimfMRIdatabase('random');

%% check each voxel tuning
%responses are missing at prior mean because we can't classify
%when direction is at prior mean. I drop those trials.
%plot switch = 1 and switch = 2
% nVox = size(d.instances,2);
% tic
% parfor voxnum = 1:nVox
%     st1{voxnum} = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch=1','myRandomCoh=0.08',d,'gridsearchMean');        
%     modes1(voxnum) = st1{voxnum}.fitvmMode; %selectivity for sw-p
%     st2{voxnum} = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch=2','myRandomCoh=0.08',d,'gridsearchMean');
%     modes2(voxnum) = st2{voxnum}.fitvmMode; %for sw-d
%     vline(225,':')
% end
% toc
% close all

[st1,modes1,st2,modes2,nVox] = slfmriGetVoxTuningForDB(d,'neuralvector',o);

%% sort voxels such that vox with selectivities closest to prior when
%switching to prior, and away from prior when switching to direction
%are first and organize all voxel plots 5 by 5 on different figures
%sorting
d2psw1 = abs(SLvectors2signedAngle(modes1,225,'polar'));
d2psw2 = -abs(SLvectors2signedAngle(modes2,225,'polar'));
[~,ix] = sortrows([d2psw1 d2psw2],[1 2]);
subploti = repmat(1:25,1,ceil(nVox/25));
subploti(nVox+1:end) = [];

for voxnum = 1:nVox   
    
    %open figure every 25 voxels
    if mod(voxnum-1,25)==0
        figure('color','w')
    end
    %plot 25 vox per figure 
    %note vox selectivity when sw-to-prior/direction
    subplot(5,5,subploti(voxnum)) 
    VoxSel1_i = modes1(ix(voxnum));
    VoxSel2_i = modes2(ix(voxnum));
    sIx = ix(voxnum); title({['Vox:' num2str(sIx)],['\color[rgb]{0 .5 1}modesP:' num2str(VoxSel1_i)],['\color[rgb]{0 0 0}modesD:' num2str(VoxSel2_i)]})%sorted index
    
    %sw-prior
    hold all
    myerrorbar(st1{sIx}.conditions(:,1),st1{sIx}.mean,'yError',...
        st1{sIx}.sem,'Symbol=o','Color=[0 .5 1]');
    plot(st1{sIx}.conditions(:,1),st1{sIx}.mean,'.','linestyle',...
        'none','color',[0 .5 1]);    
    plot(1:1:360,st1{sIx}.tuningFit,'color',[0 .5 1]); %fit
    xlim([1 360])
    
    %sw-dir
    myerrorbar(st2{sIx}.conditions(:,1),st2{sIx}.mean,'yError',...
        st2{sIx}.sem,'Symbol=o','Color=[0 0 0]');
    plot(st2{sIx}.conditions(:,1),st2{sIx}.mean,'.','linestyle',...
        'none','color',[0 0 0]);    
    plot(1:1:360,st2{sIx}.tuningFit,'color',[0 0 0]); %fit
    xlim([1 360])
    
    %legend
    if voxnum==1
        xlabel('myRandomDir')
        ylabel('Response (blue: sw-p | black: sw-d)')
    end
    %prior
    vline(225,'k:')
end



