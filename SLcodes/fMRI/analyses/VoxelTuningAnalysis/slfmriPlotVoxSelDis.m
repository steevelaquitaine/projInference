

%slfmriPlotVoxSelDis.m
%
%
% author: steeve laquitaine
%   date: 160629
%purpose: plot voxel direction selectivity distribution sorted by
%         conditions of "switching" variable (near prior vs near direction)
%
%  usage:
%
%          slfmriInitAnalysisTaskDotDirfMRI05
%          [d,o,behbySess] = slfmriGetDBoverSessions(o,varargin);
%          slfmriGetVoxTuningParams('vonMisesprior','neuralvector',d,o)
%          slfmriPlotVoxSelDis(o,VoxelsParams,'significantTuning')
%
%
% inputs :
%
%                     'all' : plot all voxels
%       'significantTuning' : plot voxels that are significantly tuned
%

function [nSignif1,nSignif2,nVox] = slfmriPlotVoxSelDis(o,VoxelsParams,voxsorting)

%number of tuned voxels
nSignif1 = [];
nSignif2 = [];
nVox = length(VoxelsParams.signifTun1);

%all voxels
if strcmp(voxsorting,'all')        
end

%keep voxels with significant tuning
if strcmp(voxsorting,'significantTuning')
    modes1 = VoxelsParams.modes1(VoxelsParams.signifTun1==1);
    modes2 = VoxelsParams.modes2(VoxelsParams.signifTun2==1);
    
    %# of significantly tuned voxels
    nSignif1 = sum(VoxelsParams.signifTun1);
    nSignif2 = sum(VoxelsParams.signifTun2);
    nVox = length(VoxelsParams.signifTun1);
end

%distribution of selectivities for switch-to-prior
figure('color','w');
subplot(121); hold all; hist(VoxelsParams.modes1,1:10:360)
slsetHistColor([.5 .5 .5],'none')
vline(o.priormean,':b'); box off
xlim([0 360]); title('Switch-to-prior')

%distribution of selectivities for switch-to-evidence
subplot(122); hold all; hist(VoxelsParams.modes2,1:10:360)
slsetHistColor([.5 .5 .5],'none')
vline(o.priormean,':b'); box off
xlim([0 360]);title('Switch-to-evidence')