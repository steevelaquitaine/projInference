

%slfMRIwrapperDisVoxSelbalDirVMfitGridS.m
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
%           %takes 32 hours (320 vox - 1000 perms)
%           slfmriInitAnalysisTaskDotDirfMRI05;
%           load d.mat
%           [voxTuning,d,o,behbySess] = slfMRIwrapperDisVoxSelbalDirVMfitGridS(o,'vmfit','significantTuning',[varargin,'db',d,'gridsearchMean']);
%
%           slfmriInitAnalysisTaskDotDirfMRI05;
%           [voxTuning,d,o,behbySess] = slfMRIwrapperDisVoxSelbalDirVMfitGridS(o,'vmfit','significantTuning',varargin);
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


function  [voxTuning,d,o,behbySess] = slfMRIwrapperDisVoxSelbalDirVMfitGridS(o,tuningType,tuningSignif,varargin)

t1 = tic;

%get the database over sessions
[d,o] = slfmriCreateOrLoadfMRIdb(o,varargin{1});

%get the distribution of voxel selectivities for equalized
%distribution of directions between "switching" conditions
[~,VoxelsParams,voxRawTunings] = slfmriGetVoxTuningParams('vonMisesprior',tuningType,d,o,varargin{1});

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
for i = 1:nVox
    voxRawTunings.st1bkp{i}.actualnVec.deg.mean = voxRawTunings.st1bkp{i}.actualnVec.datadegmean;
end
parfor voxel = 1 : nVox;
    figure(voxel)
    set(gcf,'color','w')
    slfMRIplotVoxTuningVsChance(voxRawTunings.st1bkp{voxel}.nVecs.deg.mean, ...
        voxRawTunings.st1bkp{voxel}.actualnVecStrng,...
        voxRawTunings.st1bkp{voxel}.actualnVec,voxRawTunings.st1bkp{voxel}.meantunStrngBydir,...
        voxRawTunings.st1bkp{voxel}.tunStrngCIbyDir,voxRawTunings.st1bkp{voxel}.vecdir)    
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









