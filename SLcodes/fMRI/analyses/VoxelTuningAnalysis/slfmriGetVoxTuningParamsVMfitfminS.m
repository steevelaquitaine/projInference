
%slfmriGetVoxTuningParamsVMfitfminS.m
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
%       [o,prep,VoxelsParams] = slfmriGetVoxTuningParamsVMfitfminS('myRandomDir',{'mySwitch=1','mySwitch=2'},'myRandomCoh=0.06',85,o,'useinputdb',d_dirBalbySw{coh_u==coherence},'getVoxTuning','vmfitr2cutoff=0');
%
%
%inputs:
%              var0 : tuning x space (e.g., 'myRandomDir')
%              var1 : sort tuning by variable 1 (e.g., 'mySwitch=1')
%                       var1 = {[],'mySwitch=2'}
%                       var1 = {'mySwitch=1',[]}
%                       var1 = {'mySwitch=1','mySwitch=2'}
%              var2 : another variable to select data (e.g.., 'myRandomCoh=0.06')
%          thisvar0 : get population activity for this value of var0
%
%  'loadInstancesdb': load instance database created by
%                     "slfmriGetInstancedb.m"
%       'useinputdb': use input d 
%               e.g., 'useinputdb',d
%
%     'getVoxTuning': fit voxels with von Mises to get their selectivity
%                     parameters. Time expensive.
%    'loadVoxTuning': load existing voxel selectivity parameters
%
%options:
%       'vmfitr2cutoff=0.95': set r-squared cutoff to analyse voxels which
%                             tuning was best fit with von Mises.
%               
%%outputs:
%     o.modes2u 
%     o.modes1u 
%     o.RespMeans1
%     o.RespMeans2 
%     o.RespStd2 
%     o.RespStd1 
%     o.ks1
%     o.as1
%     o.rsquared2
%     o.ks2
%     o.as2
%     o.rsquared2
%
%saves:
%       'VoxPopActivity.mat'
%       'VoxelsParams.mat'
%
%
%note: to test code with simulated instances
% d.instances = [];
% d.instances(d.mySwitch==1,:) = slSimInstancesdb(d.myRandomDir(d.mySwitch==1),repmat(225,1,360),1); %to prior
% d.instances(d.mySwitch==2,:) = slSimInstancesdb(d.myRandomDir(d.mySwitch==2),1:1:360,100);%to evidence
% save('d','d')

function [o,prep,VoxelsParams,voxTunings] = slfmriGetVoxTuningParamsVMfitfminS(var0,var1,var2,thisvar0,o,varargin)

%preprocess
prep = getPreprocessing(varargin{:});

%load instance database
cd(o.dataPath)
d = slfmriLoadInstancedb(varargin{:});

%--------------- Voxels tuning -----------------
%when we need to get voxels tuning (takes some time)
%more for many voxels
if any(strcmp(varargin,'getVoxTuning'))
    fprintf('%s \n','(slfmriGetVoxTuningParamsVMfitfminS) Getting voxel tunings.')
    nVox = size(d.instances,2);
    
    %------- variable 1 value 1 ----------
    %e.g., switch-to-prior or coherence = 0.06
    %best von mises tuning fit for each voxel    
    %check enough sample in var0 for fitting circular to data
    if length(unique(d.(var0)))<3
        fprintf('%s \n','(slfmriGetVoxTuningParamsVMfitfminS) Stopped. Not enough sample in x space (n<3) to fit von Mises to data. Fitting will produce meaningless results.')
        keyboard        
    end
    
    %fit
    modes1 = nan(nVox,1); ks1 = nan(nVox,1); as1 = nan(nVox,1); rsquared1 = nan(nVox,1);
    if ~isempty(var1{1})
        for voxnum = 1 : nVox
            fprintf('-')
            st1 = slfmrigetVoxTuning(voxnum,var0,var1{1},var2,d,'fminsearch');
            modes1(voxnum) = st1.fitvmMode;
            ks1(voxnum)    = st1.fitvmK;
            as1(voxnum)    = st1.peakAmp;
            signifTun1(voxnum) = st1.signifTun;                        
            rsquared1(voxnum) = st1.rsquared;
            st1bkp{voxnum} = st1;            
        end
        figure; plot(rsquared1)
    end
    
    %------- variable 1 value 2 ----------
    %e.g., switch-to prior or coherence = 0.12        
    %best von mises fit of each voxel tuning
    modes2 = nan(nVox,1); ks2= nan(nVox,1); as2 = nan(nVox,1); rsquared2 = nan(nVox,1);    
    if ~isempty(var1{2})
        parfor  voxnum = 1 : nVox
            fprintf('-')
            st2 = slfmrigetVoxTuning(voxnum,var0,var1{2},var2,d,'fminsearch');
            modes2(voxnum) = st2.fitvmMode;
            ks2(voxnum)    = st2.fitvmK;
            as2(voxnum)    = st2.peakAmp;
            signifTun2(voxnum) = st2.signifTun;                        
            rsquared2(voxnum) = st2.rsquared;
            st2bkp{voxnum} = st2;
        end
        figure; plot(rsquared2)        
    end
    
    %save params
    voxTunings.st1bkp = st1bkp;
    voxTunings.st2bkp = st2bkp;
    save('voxTunings','voxTunings')
    
    %save params
    VoxelsParams.modes1 = modes1;
    VoxelsParams.tunStrng1 = ks1;
    VoxelsParams.as1 = as1;
    VoxelsParams.signifTun1 = signifTun1;
    VoxelsParams.rsquared1 = rsquared1;
    VoxelsParams.modes2 = modes2;
    VoxelsParams.tunStrng2 = ks2;
    VoxelsParams.as2 = as2;
    VoxelsParams.signifTun2 = signifTun2;
    VoxelsParams.rsquared2 = rsquared2;
    VoxelsParams.description = {'1/2: switch-to-prior/dir, modes: voxels selectivities'};
        
    %save('VoxelsParams','modes1','as1','signifTun1','modes2','as2','signifTun2','description')    
    save('VoxelsParams','VoxelsParams')    
        
    %when vox tuning data already exist (faster)
elseif any(strcmp(varargin,'loadVoxTuning'))
    load VoxelsParams
end
    

%-----------------Nested ----------
%get preprocessing analysis
function prep = getPreprocessing(varargin)

%get analyses
vmfitr2cutoff = [];
getArgs(varargin,{'vmfitr2cutoff=[]'});

%store
prep.vmfitr2cutoff = vmfitr2cutoff;


