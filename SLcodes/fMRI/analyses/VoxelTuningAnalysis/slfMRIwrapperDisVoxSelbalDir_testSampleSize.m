
%slfMRIwrapperDisVoxSelbalDir_testSampleSize.m
%author: steeve laquitaine
%  date: 160423


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


%% Equalize distribution of directions across switching variable
du = unique(dcut.myRandomDir);
nDu = length(du);

%coherence 6%: get trials to keep to balance direction by switch variable
ixDiriBal_sw1 = [];
ixDiriBal_sw2 = [];
parfor i = 1 : nDu
    %trials this direction & condition
    ixDiri_Sw1 = find(dcut.mySwitch==1 & dcut.myRandomCoh==0.06 & dcut.myRandomDir==du(i));
    ixDiri_Sw2 = find(dcut.mySwitch==2 & dcut.myRandomCoh==0.06 & dcut.myRandomDir==du(i));
    %equalize number of occurence for each direction
    %by removing additional directions
    minNb = min([length(ixDiri_Sw1); length(ixDiri_Sw2)]);    
    %we keep the first n min trials
    ixDiriBal_sw1 = [ixDiriBal_sw1; ixDiri_Sw1(1:minNb)];
    ixDiriBal_sw2 = [ixDiriBal_sw2; ixDiri_Sw2(1:minNb)];    
end
%re-build the database
d_dirBalbySw006 = dcut;
fieldnm = fieldnames(dcut);
ixAll = sort([ixDiriBal_sw1;ixDiriBal_sw2]);
for i = 1:length(fieldnames(dcut))
    d_dirBalbySw006.(fieldnm{i}) = dcut.(fieldnm{i})(ixAll,:);
end

%coherence 12%: get trials to keep to balance direction by switch variable
ixDiriBal_sw1 = [];
ixDiriBal_sw2 = [];
parfor i = 1 : nDu
    ixDiri_Sw1 = find(dcut.mySwitch==1 & dcut.myRandomCoh==0.12 & dcut.myRandomDir==du(i));
    ixDiri_Sw2 = find(dcut.mySwitch==2 & dcut.myRandomCoh==0.12 & dcut.myRandomDir==du(i));   
    minNb = min([length(ixDiri_Sw1); length(ixDiri_Sw2)]);    
    ixDiriBal_sw1 = [ixDiriBal_sw1; ixDiri_Sw1(1:minNb)];
    ixDiriBal_sw2 = [ixDiriBal_sw2; ixDiri_Sw2(1:minNb)];    
end
d_dirBalbySw012 = dcut;
fieldnm = fieldnames(dcut);
ixAll = sort([ixDiriBal_sw1;ixDiriBal_sw2]);
for i = 1:length(fieldnames(dcut))
    d_dirBalbySw012.(fieldnm{i}) = dcut.(fieldnm{i})(ixAll,:);
end


%% plot distribution of directions sorted by switching variable
figure('color','w'); 
%6%
subplot(1,2,1); hold all
h1 = hist(d_dirBalbySw006.myRandomDir(d_dirBalbySw006.mySwitch==1),1:1:360);
bar([1:1:360]-4,h1,10,'facecolor',[0 0.5 .7],'lineStyle','none')
h2 = hist(d_dirBalbySw006.myRandomDir(d_dirBalbySw006.mySwitch==2),1:1:360);
bar([1:1:360]+4,h2,10,'facecolor','k','lineStyle','none')
ylabel('Count'); xlabel('Motion direction (deg)'); xlim([0 360])
legend('Sw-p','Sw-d'); alpha(0.9)
set(gca,'xtick',unique(d_dirBalbySw006.myRandomDir),'xticklabel',unique(d_dirBalbySw006.myRandomDir))
vline(225,':k') ; title('Direction dist. by switching for 6% coh'); box off
%12%
subplot(1,2,2); hold all
h3 = hist(d_dirBalbySw012.myRandomDir(d_dirBalbySw012.mySwitch==1),1:1:360);
bar([1:1:360]-4,h3,10,'facecolor',[0 0.5 .7],'lineStyle','none')
h4 = hist(d_dirBalbySw012.myRandomDir(d_dirBalbySw012.mySwitch==2),1:1:360);
bar([1:1:360]+4,h4,10,'facecolor','k','lineStyle','none')
ylabel('Count'); xlabel('Motion direction (deg)'); xlim([0 360]); legend('Sw-p','Sw-d')
alpha(0.9); set(gca,'xtick',unique(d_dirBalbySw012.myRandomDir),'xticklabel',unique(d_dirBalbySw012.myRandomDir))
vline(225,':k') ;title('For 12% coh') ;box off




%% distribution of all voxel selectivities
o = slfmriGetVoxPopActivity2('myRandomDir',{'mySwitch=1','mySwitch=2'},'myRandomCoh=0.06',85,o,'useinputdb',d_dirBalbySw006,'getVoxTuning','vmfitr2cutoff=0');

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




