
%slfMRIwrapperPlotVoxTuningDB_testSampleSize.m
%
%
% author: steeve laquitaine
%   date: 160310
%purpose: plot database of voxel tuning (e.g., to stimulus motion direction)
%
%Description: 
%   1) initialize data and analysis info
%   2) make database (take half the sample)
%   3) vizualize sorted voxel tunings by switch variable

%initialize analysis and make fMRI database
slfmriInitAnalysisTaskDotDirfMRI05
[d,o] = slfMRImakeDatabase(o,varargin);

%% cut # of trials by half
%backup d
dcut = d;

%find trials each direction
ix15_06 = find(d.myRandomDir==15 &d.myRandomCoh==0.06);
ix85_06 = find(d.myRandomDir==85 &d.myRandomCoh==0.06);
ix155_06 = find(d.myRandomDir==155 &d.myRandomCoh==0.06);
ix225_06 = find(d.myRandomDir==225 &d.myRandomCoh==0.06);
ix295_06 = find(d.myRandomDir==295 &d.myRandomCoh==0.06);

%reduce to half the trials
ix15_06half = ix15_06(1:floor(length(ix15_06)/2));
ix85_06half = ix85_06(1:floor(length(ix85_06)/2));
ix155_06half = ix155_06(1:floor(length(ix155_06)/2));
ix225_06half = ix225_06(1:floor(length(ix225_06)/2));
ix295_06half = ix295_06(1:floor(length(ix295_06)/2));
ixhalf = [ix15_06half;ix85_06half;ix155_06half;ix225_06half;ix295_06half];

%update the database
dcut.instances = dcut.instances(ixhalf,:);
dcut.stimvols = dcut.stimvols(ixhalf,:);
dcut.myRandomCoh = dcut.myRandomCoh(ixhalf,:);
dcut.myRandomDir = dcut.myRandomDir(ixhalf,:);
dcut.mySwitch = dcut.mySwitch(ixhalf,:);


%% check the tuning fit of individual voxel
%responses are missing at prior mean because it is impossible to classify
%estimates as switching to prior or direction when motion direction is
%displayed at prior mean. Thus those trials are removed.
%plot switch = 1 and switch = 2
nVox = size(dcut.instances,2);
parfor voxnum = 1:nVox
    st1{voxnum} = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch=1','myRandomCoh=0.06',dcut,'fminsearch');        
    st2{voxnum} = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch=2','myRandomCoh=0.06',dcut,'fminsearch');
    vline(225,':')
    modes1(voxnum) = st1{voxnum}.fitvmMode; %selectivity for sw-p
    modes2(voxnum) = st2{voxnum}.fitvmMode; %for sw-d
end
close all

%% sort voxels such that vox with selectivities closest to prior when
%switching to prior and furthest from prior when switching to direction
%are first and plot
%plot the many voxels 5 by 5 on many figures
figure('color','w')
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
    subplot(5,5,subploti(voxnum))    
    sIx = ix(voxnum); title(['Vox:' num2str(sIx)])%sorted index
    
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




