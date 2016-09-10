

%slfMRIwrapperDisVoxSelbalDir.m
%
%
% author: steeve laquitaine
%   date: 160314
%purpose: plot distribution of voxel selectivities sorted by when subject 
%         switches to prior or evidence for equalized motion direction 
%         distributions between switching conditions
%
%   usage :
% 
%           slfmriInitAnalysisTaskDotDirfMRI05;
%           [voxTuning,d,o,behbySess] = slfMRIwrapperDisVoxSelbalDir(o,'neuralvector','significantTuning',varargin);
%
%           slfmriInitAnalysisTaskDotDirfMRI05;
%           [voxTuning,d,o,behbySess] = slfMRIwrapperDisVoxSelbalDir(o,'vmfit','significantTuning',varargin);
%
%
% inputs :
%
%       'rawpeak': selectivity is the direction that produces the maximal response
%         'vmfit': selectivity is the mode of a fitted von Mises
%  'neuralvector': selectivity is the voxel direction vector 
%                  (voxel responses weight their associated directions Georgopoulos)
%
%
%
%           'significantTuning'
%
%note : tried with 1 millions permutations: the computer crashes. Warning 
%       that it can't save voxRawTunings (probably too heavy)


function  [voxTuning,d,o,behbySess] = slfMRIwrapperDisVoxSelbalDir(o,tuningType,tuningSignif,varargin)

t1 = tic;

%get the database over sessions
[d,o,behbySess] = slfmriGetDBoverSessions(o,varargin);

%get the distribution of voxel selectivities for equalized
%distribution of directions between "switching" conditions
[~,VoxelsParams,voxRawTunings] = slfmriGetVoxTuningParams('vonMisesprior',tuningType,d,o);

%plot voxel selectivity distributions
[nSignif1,nSignif2,nVox] = slfmriPlotVoxSelDis(o,VoxelsParams,tuningSignif);
subplot(1,2,1);
title(['Switch-to-prior (' num2str(nSignif1) '/' num2str(nVox) ' signif tuned voxels)'])
xlabel({'Voxel selectivities',['(' tuningType ', ' tuningSignif ',' o.myROIname{1} ')' ]})
subplot(1,2,2);
title(['Switch-to-evidence (' num2str(nSignif2) '/' num2str(nVox) ' signif tuned voxels)'])
xlabel({'Voxel selectivities',['(' tuningType ', ' tuningSignif ',' o.myROIname{1} ')' ]})

%plot each voxel observed tuning vs. its null tunings distribution from
%permuted responses
parfor voxel = 1 : nVox;
    figure(voxel)
    set(gcf,'color','w')
    slfMRIplotVoxTuningVsChance(voxRawTunings.st1bkp{voxel}.nVecs.deg.mean, ...
        voxRawTunings.st1bkp{voxel}.actualnVeclen,...
        voxRawTunings.st1bkp{voxel}.actualnVec,voxRawTunings.st1bkp{voxel}.meanVlenBydir,...
        voxRawTunings.st1bkp{voxel}.VlenCIbyDir,voxRawTunings.st1bkp{voxel}.vecdir)    
    slgraphicsSaveEPS
    close;
end
toc(t1)

%outputs
voxTuning.VoxelsParams = VoxelsParams;
voxTuning.voxRawTunings = voxRawTunings;
voxTuning.nSignif1 = nSignif1;
voxTuning.nSignif2 = nSignif2;
voxTuning.nVox = nVox;
o.duration = toc(t1);
fprintf('%s %i %s \n','Took',o.duration,'sec')









