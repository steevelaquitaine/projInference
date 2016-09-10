
%slfmriGetVoxPopActivity2.m
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
%       o = slfmriGetVoxPopActivity2('myRandomDir',{'mySwitch=1','mySwitch=2'},'myRandomCoh=0.06',85,o,'useinputdb',d_dirBalbySw{coh_u==coherence},'getVoxTuning','vmfitr2cutoff=0');
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

function [o,prep] = slfmriGetVoxPopActivity2(var0,var1,var2,thisvar0,o,varargin)

%preprocess
prep = getPreprocessing(varargin{:});

%load instance database with variables
cd(o.dataPath)

%load instance database created by slfmriGetInstancedb
if any(strcmp(varargin,'loadInstancesdb'))
    d = importdata('d.mat');
elseif any(strcmp(varargin,'useinputdb'))
    fprintf('%s \n','(slfmriGetVoxPopActivity2) Using input database.')
    d = varargin{find(strcmp(varargin,'useinputdb'))+1};
end

%--------------- Voxels tuning -----------------
%when we need to get voxels tuning (takes some time)
%more for many voxels
if any(strcmp(varargin,'getVoxTuning'))
    fprintf('%s \n','(slfmriGetVoxPopActivity2) Getting voxel tunings.')
    nVox = size(d.instances,2);
    
    %------- variable 1 value 1 ----------
    %e.g., switch-to-prior or coherence = 0.06
    %preallocate for speed
    modes1 = nan(nVox,1);
    ks1 = nan(nVox,1);
    as1 = nan(nVox,1);
    rsquared1 = nan(nVox,1);
    
    %best von mises tuning fit for each voxel    
    %check enough sample in var0 for fitting circular to data
    if length(unique(d.(var0)))<3
        fprintf('%s \n','(slfmrigetVoxTuning) Stopped. Not enough sample in x space (n<3) to fit von Mises to data. Fitting will produce meaningless results.')
        keyboard        
    end
    
    %fit
    if ~isempty(var1{1})
        parfor voxnum = 1 : nVox
            fprintf('-')
            st1 = slfmrigetVoxTuning(voxnum,var0,var1{1},var2,d,'gridsearchMean');
            %backup params
            modes1(voxnum) = st1.fitvmMode;
            ks1(voxnum) = st1.fitvmK;
            as1(voxnum) = st1.peakAmp;
            %check goodness of fit
            rsquared1(voxnum) = st1.rsquared;
        end
        figure; plot(rsquared1)
    end
    %------- variable 1 value 2 ----------
    %e.g., switch-to prior or coherence = 0.12
    %preallocate for speed
    modes2 = nan(nVox,1);
    ks2    = nan(nVox,1);
    as2    = nan(nVox,1);
    rsquared2 = nan(nVox,1);
        
    %best von mises fit of each voxel tuning
    if ~isempty(var1{2})
        parfor  voxnum = 1 : nVox
            fprintf('-')
            %get tunings
            %st2 = slfmrigetVoxTuning(voxnum,var0,var1{2},d,'nofit');
            st2 = slfmrigetVoxTuning(voxnum,var0,var1{2},var2,d,'gridsearchMean');
            %backup params
            modes2(voxnum) = st2.fitvmMode;
            ks2(voxnum)    = st2.fitvmK;
            as2(voxnum)    = st2.peakAmp;
            %check goodness of fit
            rsquared2(voxnum) = st2.rsquared;
        end
        figure; plot(rsquared2)        
    end
    
    %save params
    description = {'1/2: switch-to-prior/dir, modes: voxels selectivities',...
        ;'ks: best fit von mises concentration parameter'};
    save('VoxelsParams','modes1','ks1','as1','rsquared1',...
        'modes2','ks2','as2','rsquared2','description')    
    
    %when vox tuning data already exist (faster)
elseif any(strcmp(varargin,'loadVoxTuning'))
    load VoxelsParams
end

%--------------- Population activity by var1 -----------------
%check if any voxel at all pass chosen r-squared cutoff
% nHighr2 = sum(rsquared1>=prep.vmfitr2cutoff);
% if nHighr2 == 0
%     fprintf('%s \n','No voxel passed the r-squared cutoff. Lower it down.')
%     return
% end
% nHighr2 = sum(rsquared2>=prep.vmfitr2cutoff);
% if nHighr2 == 0
%     fprintf('%s \n','No voxel passed the r-squared cutoff. Lower it down.')
%     return
% end

%get responses by selectivity for switch-to-prior
figure('color','w')

%unique selectivities when switch-to-prior
modes1u = unique(modes1);
nmodes1u = length(modes1u);
thisCond =[];
% thisCond = d.myRandomDir==thisvar0 & d.mySwitch==1;
thisCond = d.myRandomDir==thisvar0 & d.mySwitch==2;

%preallocate for speed
thisMode = [];
thisModeRes1 = cell(nmodes1u,1); 
RespMeans1   = nan(nmodes1u,1); 
RespStd1     = nan(nmodes1u,1); 

%for each selectivity in tuning x space
for i = 1 : nmodes1u
    
    %voxel selectivities for voxels tuning that
    %satisfy r2cutoff (e.r., high r2 indicates good
    %von Mises fit of voxel tuning)    
    if ~isempty(prep.vmfitr2cutoff)
        thisMode = modes1==modes1u(i) & rsquared1>=prep.vmfitr2cutoff;
    else        
        thisMode = modes1==modes1u(i);
    end    
    
    %responses (averaged of repeats)
    if any(thisCond==1)
        if any(thisMode==1)
            thisModeRes1{i} = mean(d.instances(thisCond,thisMode),1);
            %%plot each voxel
            %hold on; plot(modes1u(i),thisModeRes1{i},'.');            
        else
            fprintf('%s \n','No condition found')
        end
    end
    %mean response over vox with same selectivities and std 
    RespMeans1(i) = mean(thisModeRes1{i});
    RespStd1(i) = std(thisModeRes1{i});  
    
    %save plot stats
    RespnumVox1(i) = length(thisModeRes1{i});
end
%plot population activity by selectivity for value 1 of var1
title(var1{1})
hold on; myerrorbar(modes1u,RespMeans1,'yError',RespStd1,'Color=[0 .5 .9]')
xlim([0 360])
vline(thisvar0,'k--')
vline(o.priormean,'b--')


%response by voxel selectivities sorted by var1 and var0
%unique selectivities when "switch-to-direction"
%selectivities for value 1 of var 1 
modes2u  = unique(modes2);
nmodes2u = length(modes2u);
figure('color','w')

%set this value of var0 and value 2 of var1
thisCond = d.myRandomDir==thisvar0 & d.mySwitch==2;

%preallocate for speed
thisMode = [];
thisModeRes2 = cell(nmodes2u,1); 
RespMeans2   = nan(nmodes2u,1); 
RespStd2     = nan(nmodes2u,1); 

%for each voxel selectivity found in tuning x space
for i = 1 : nmodes2u
    
    %this selectivity for voxels tuning that
    %satisfy r2cutoff (e.r., high r2 indicates good
    %von Mises fit of voxel tuning)
    if ~isempty(prep.vmfitr2cutoff)
        thisMode = modes2==modes2u(i) & rsquared2>=prep.vmfitr2cutoff;
    else        
        thisMode = modes2==modes2u(i);
    end
    
    %check condition exists
    if any(thisCond==1)
        if any(thisMode==1)
            
            %responses (averaged over repeats)
            thisModeRes2{i} = mean(d.instances(thisCond,thisMode),1);
            
        else
            fprintf('%s \n','No condition found')
        end
    end
    %mean response and std by selectivity (averaged over voxels)
    RespMeans2(i) = mean(thisModeRes2{i});
    RespStd2(i) = std(thisModeRes2{i}); 
    
    %save plot stats
    RespnumVox2(i) = length(thisModeRes2{i});
end
%plot population activity by selectivity for value 2 of var1
title(var1{2})
hold on; myerrorbar(modes2u,RespMeans2,'yError',RespStd2)
xlim([0 360])
vline(thisvar0,'k--')
vline(o.priormean,'b--')

%outputs
o.modes2u    = modes2u;
o.modes1u    = modes1u;
o.RespMeans1 = RespMeans1;
o.RespMeans2 = RespMeans2;
o.RespStd2   = RespStd2;
o.RespStd1   = RespStd1;
o.RespnumVox1 = RespnumVox1;
o.RespnumVox2 = RespnumVox2;
o.ks1         = ks1;
o.as1         = as1;
o.rsquared1   = rsquared1;
o.ks2         = ks2;
o.as2         = as2;
o.rsquared2   = rsquared2;

%describe output variables
o.description = {'modes2u','modes1u','RespMeans1','RespMeans2',...
    'RespStd2','RespStd1',...
    'RespnumVox1 : # of voxels satisfying conditions by selectivity when variable 1 val 1',...
    'RespnumVox1 : # of voxels satisfying conditions by selectivity when variable 1 val 2'};

%save results
save('VoxPopActivity','modes2u','modes1u','RespMeans1','RespMeans2',...
    'RespnumVox1','RespnumVox2','RespStd2','RespStd1',...
    'rsquared1','rsquared2','prep','o')

    

%-----------------Nested ----------
%get preprocessing analysis
function prep = getPreprocessing(varargin)

%get analyses
vmfitr2cutoff = [];
getArgs(varargin,{'vmfitr2cutoff=[]'});

%store
prep.vmfitr2cutoff = vmfitr2cutoff;


