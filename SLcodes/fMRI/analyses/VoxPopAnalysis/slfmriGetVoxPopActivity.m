
%slfmriGetVoxPopActivity.m
%
%
% author: steeve laquitaine
%   date: 160117
% pupose: get voxel population (response by selectivities) sorted by a
%         variable (for now only switching)
%
% usage:
%
%       o.dataPath = '/Volumes/Transcend/data/sltaskdotdirfmri05/slStckAnalyses/MotionComp/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/IPS/';
%       [o,prep] = slfmriGetVoxPopActivity('myRandomDir',{'mySwitch=1','mySwitch=2'},15,o,'loadVoxTuning','loadInstancesdb')
%
%inputs:
%              var0 : tuning x space (e.g., 'myRandomDir')
%              var1 : sort tuning by variable 1 (e.g., 'mySwitch=1')
%          thisvar0 : get population activity for this value of var0
%
%  'loadInstancesdb': load instance database created by
%                     "slfmriGetInstancedb.m"
%     'getVoxTuning': fit voxels with von Mises to get their selectivity
%                     parameters. Time expensive.
%    'loadVoxTuning': load existing voxel selectivity parameters
%
%options:
%       'vmfitr2cutoff=0.95': set r-squared cutoff to analyse voxels which
%                             tuning was best fit with von Mises.
%               
%%outputs
%     o.modesDu 
%     o.modesPu 
%     o.RespMeansP 
%     o.RespMeansD 
%     o.RespStdD 
%     o.RespStdP 
%     o.ksP
%     o.asP
%     o.rsquaredP 
%     o.ksD
%     o.asD
%     o.rsquaredD 
%
function [o,prep] = slfmriGetVoxPopActivity(var0,var1,thisvar0,o,varargin)

%preprocess
prep = getPreprocessing(varargin{:});

%load instance database with variables
cd(o.dataPath)

%load instance database created by slfmriGetInstancedb
if any(strcmp(varargin,'loadInstancesdb'))
    d = importdata('d.mat');
end

%when we need to get voxels tuning (takes some time)
%more for many voxels
if any(strcmp(varargin,'getVoxTuning'))
    nVox = size(d.instances,2);
    %------- switch-to-prior ----------
    %preallocate for speed
    modesP = nan(nVox,1);
    ksP = nan(nVox,1);
    asP = nan(nVox,1);
    rsquaredP = nan(nVox,1);
    %best von mises tuning fit for each voxel
    parfor  voxnum = 1 : nVox        
        stP = slfmrigetVoxTuning(voxnum,var0,var1{1},d,'fminsearch');
        %backup params
        modesP(voxnum) = stP.fitvmMode;
        ksP(voxnum) = stP.fitvmK;
        asP(voxnum) = stP.peakAmp;
        %check goodness of fit
        rsquaredP(voxnum) = stP.rsquared;
    end
    figure; plot(rsquaredP)

    %------- switch-to-direction ----------
    %preallocate for speed
    modesD = nan(nVox,1);
    ksD = nan(nVox,1);
    asD = nan(nVox,1);
    rsquaredD = nan(nVox,1);
    %best von mises fit of each voxel tuning
    parfor  voxnum = 1 : nVox
        stD = slfmrigetVoxTuning(voxnum,var0,var1{2},d,'fminsearch');
        %backup params
        modesD(voxnum) = stD.fitvmMode;
        ksD(voxnum) = stD.fitvmK;
        asD(voxnum) = stD.peakAmp;
        %check goodness of fit
        rsquaredD(voxnum) = stD.rsquared;
    end
    figure; plot(rsquaredD)
    
    %save params
    save('VoxelsParams','modesP','ksP','asP','rsquaredP',...
        'modesD','ksD','asD','rsquaredP')    
    
    %when vox tuning data already exist (faster)
elseif any(strcmp(varargin,'loadVoxTuning'))
    load VoxelsParams
end

%--------------- Population activity by var1 -----------------
%response by voxel selectivities sorted by var1 and var0
%unique selectivities when "switch-to-direction"
%selectivities for value 1 of var 1 
modesDu = unique(modesD);
nmodesDu = length(modesDu);
%get responses by selectivity for value 1 of var 1 
figure('color','w')
%set this value of var0 and value 2 of var1
sw = 2;
thisCond = d.myRandomDir==thisvar0 & d.mySwitch==sw;
%preallocate for speed
thisMode = [];
thisModeResD = cell(nmodesDu,1); 
RespMeansD = nan(nmodesDu,1); 
RespStdD = nan(nmodesDu,1); 

%for each voxel selectivity found in tuning x space
for i = 1 : nmodesDu
    %this selectivity for voxels tuning that
    %satisfy r2cutoff (e.r., high r2 indicates good
    %von Mises fit of voxel tuning)
    thisMode = modesD==modesDu(i) & rsquaredD>=prep.vmfitr2cutoff;
    %responses (averaged over repeats)
    if any(thisCond==1)
        if any(thisMode==1)
            thisModeResD{i} = mean(d.instances(thisCond,thisMode),1);
            %%plot each voxel
            %hold on; plot(modesDu(i),thisModeResD{i},'.');
        else
            fprintf('%s \n','No condition found')
        end
    end
    %mean response and std by selectivity
    RespMeansD(i) = mean(thisModeResD{i});
    RespStdD(i) = std(thisModeResD{i}); 
    
    %save plot stats
    RespnumVoxD(i) = length(thisModeResD{i});
end
%plot population activity by selectivity for value 2 of var1
hold on; myerrorbar(modesDu,RespMeansD,'yError',RespStdD)
xlim([0 360])
vline(thisvar0,'k--')
vline(225,'b--')

%check if any voxel at all pass chosen r-squared cutoff
nHighr2 = sum(rsquaredP>=prep.vmfitr2cutoff);
if nHighr2 == 0
    fprintf('%s \n','No voxel passed the r-squared cutoff. Lower it down.')
    return
end
nHighr2 = sum(rsquaredD>=prep.vmfitr2cutoff);
if nHighr2 == 0
    fprintf('%s \n','No voxel passed the r-squared cutoff. Lower it down.')
    return
end

%get responses by selectivity for switch-to-prior
figure('color','w')
%set value 1 of var1
sw = 1;
%unique selectivities when switch-to-prior
modesPu = unique(modesP);
nmodesPu = length(modesPu);
thisCond =[];
thisCond = d.myRandomDir==thisvar0 & d.mySwitch==sw;
%preallocate for speed
thisMode = [];
thisModeResP = cell(nmodesPu,1); 
RespMeansP = nan(nmodesPu,1); 
RespStdP = nan(nmodesPu,1); 
%for each selectivity in tuning x space
for i = 1 : nmodesPu
    %voxel selectivities for voxels tuning that
    %satisfy r2cutoff (e.r., high r2 indicates good
    %von Mises fit of voxel tuning)
    thisMode = modesP==modesPu(i) & rsquaredP>=prep.vmfitr2cutoff;    
    %responses (averaged of repeats)
    if any(thisCond==1)
        if any(thisMode==1)
            thisModeResP{i} = mean(d.instances(thisCond,thisMode),1);
            %%plot each voxel
            %hold on; plot(modesPu(i),thisModeResP{i},'.');            
        else
            fprintf('%s \n','No condition found')
        end
    end
    %mean response and std by selectivity
    RespMeansP(i) = mean(thisModeResP{i});
    RespStdP(i) = std(thisModeResP{i});  
    
    %save plot stats
    RespnumVoxP(i) = length(thisModeResP{i});
end
%plot population activity by selectivity for value 1 of var1
hold on; myerrorbar(modesPu,RespMeansP,'yError',RespStdP,'Color=[0 .5 .9]')
xlim([0 360])
vline(thisvar0,'k--')
vline(225,'b--')

%outputs
o.modesDu    = modesDu;
o.modesPu    = modesPu;
o.RespMeansP = RespMeansP;
o.RespMeansD = RespMeansD;
o.RespStdD   = RespStdD;
o.RespStdP   = RespStdP;
o.RespnumVoxP = RespnumVoxP;
o.RespnumVoxD = RespnumVoxD;
o.ksP         = ksP;
o.asP         = asP;
o.rsquaredP   = rsquaredP;
o.ksD         = ksD;
o.asD         = asD;
o.rsquaredD   = rsquaredD;

%describe output variables
o.description = {'modesDu','modesPu','RespMeansP','RespMeansD',...
    'RespStdD','RespStdP',...
    'RespnumVoxP : # of voxels satisfying conditions by selectivity when switch to prior',...
    'RespnumVoxP : # of voxels satisfying conditions by selectivity when switch to direction'};

%save results
save('VoxPopActivity','modesDu','modesPu','RespMeansP','RespMeansD',...
    'RespnumVoxP','RespnumVoxD','RespStdD','RespStdP',...
    'rsquaredP','rsquaredD','prep','o')

    

%-----------------Nested ----------
%get preprocessing analysis
function prep = getPreprocessing(varargin)

%get analyses
getArgs(varargin,{'vmfitr2cutoff=0'});

%store
prep.vmfitr2cutoff = vmfitr2cutoff;


