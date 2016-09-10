
%slfmriGetVoxTuningParams.m
%
% author: steeve laquitaine
%   date: 160629
%purpose: wrapper to get voxel direction tuning parameters quantified by
%         - calculating voxel tuning by taking the argmax response
%         - fitting a von Mises function to voxels responses or
%         - calculating voxel directional vector
%
%usage:
%
%       o = slfmriGetVoxTuningParams('vonMisesprior','rawpeak',d,o)
%
%
% Inputs :
%
% 'vonMisesprior': session with von Mises distribution of motion
%                       direction
%  'uniformPrior': session with uniform distribution of motion
%                       directions
% 
%       'rawpeak': selectivity is the direction that produces the maximal response
%         'vmfit': selectivity is the mode of a fitted von Mises
%  'neuralvector': selectivity is the voxel direction vector 
%                  (voxel responses weight their associated directions Georgopoulos)
%
%              d : database of voxel response instances and associated
%                  coherence, direction, switching and stimvols variables
%                  retrieved by "slfmriGetDBoverSessions.m"
%
%              o : analysis parameters set by "slfmriInitAnalysisTaskDotDirfMRI05"
%
%   VoxelsParams : voxel parameters: selectivities etc...


function [o,VoxelsParams,voxRawTunings] = slfmriGetVoxTuningParams(priorType,tuningType,d,o,varargin)

VoxelsParams = [];

%case voxels tuning is fit with von Mises
if strcmp(tuningType,'vmfit')
    if strcmp(priorType,'vonMisesprior');        
        [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist00(d,o);
        varargin{1} = {'useinputdb',d_dirBalbySw{coh_u==o.coh2plot},'getVoxTuning','vmfitr2cutoff=0',varargin{1}{:}};
        [o,~,VoxelsParams,voxRawTunings] = slfmriGetVoxTuningParamsVMfit('myRandomDir',{'mySwitch=1','mySwitch=2'},['myRandomCoh=' num2str(o.coh2plot)],o.motDir,o,varargin{1});        
    end    
    if strcmp(priorType,'uniformPrior');
        [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist(d,o);
        o = slfmriGetVoxTuningParamsVMfit('myRandomDir',{[],'mySwitch=2'},'myRandomCoh=1',85,o,'useinputdb',d_dirBalbySw{coh_u==o.coh2plot},'getVoxTuning','vmfitr2cutoff=0');
    end
end

%case voxels tuning is just the direction that elicits the 
%largest voxel response
if strcmp(tuningType,'rawpeak')
    if strcmp(priorType,'vonMisesprior');
        [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist00(d,o);
        o = slfmriGetVoxTuningRawPeak('myRandomDir',{'mySwitch=1','mySwitch=2'},['myRandomCoh=' num2str(o.coh2plot)],o.motDir,o,'useinputdb',d_dirBalbySw{coh_u==o.coh2plot},'getVoxTuning','vmfitr2cutoff=0');
    end    
    if strcmp(priorType,'uniformPrior');
        [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist(d,o);
        o = slfmriGetVoxTuningRawPeak('myRandomDir',{[],'mySwitch=2'},'myRandomCoh=1',85,o,'useinputdb',d_dirBalbySw{coh_u==o.coh2plot},'getVoxTuning','vmfitr2cutoff=0');
    end
end

%case voxels tuning is just neural vector
if strcmp(tuningType,'neuralvector')
    if strcmp(priorType,'vonMisesprior');
        [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist00(d,o);
        [o,~,VoxelsParams,voxRawTunings] = slfmriGetVoxTuningNeuralVector('myRandomDir',{'mySwitch=1','mySwitch=2'},['myRandomCoh=' num2str(o.coh2plot)],o.motDir,o,'useinputdb',d_dirBalbySw{coh_u==o.coh2plot},'getVoxTuning','vmfitr2cutoff=0');
    end    %save('VoxelsParams','modes1','as1','signifTun1','modes2','as2','signifTun2','description')    
    if strcmp(priorType,'uniformPrior');
        [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist(d,o);
        [o,~,VoxelsParams,voxRawTunings] = slfmriGetVoxTuningNeuralVector('myRandomDir',{[],'mySwitch=2'},'myRandomCoh=1',85,o,'useinputdb',d_dirBalbySw{coh_u==o.coh2plot},'getVoxTuning','vmfitr2cutoff=0');
    end
end