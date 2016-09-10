
%SLsimPopVoxAnalforBayesianIntegVsCompetition.m
%
% author: steeve laquitaine
%   date: 141010
%purpose: predictions about change in BOLD in Bayesian posterior 
%         encodding-area


function SLsimPopVoxAnalforBayesianIntegVsCompetition

%--------------------
%Bayesian integration
%--------------------
close all
figure(1)

%voxels tuning
scale2unit = 10;
tuningVox135 = scale2unit*vmPdfs(1:1:360,135,2.4,[]);
tuningVox315 = scale2unit*vmPdfs(1:1:360,315,2.4,[]);

%plot tuning
drawVoxTuning(tuningVox135,tuningVox315,[693 756 617 324],'linestyle','-','color',[0.5 .1 0])
title({'Bayesian integration: stronger priors', 'induce larger shifts in voxels seletivities matching',...
    ' perceptual biases'})

%plot qualitative simulation of prior-induced shift
%for simplicity I model the shift by applying a gaussian prior gain to
%voxels activities. Multiplying by a gaussian prior gain produces narrower
%tunings for larger prior gains is "as if" we assumed homeostasis that is 
%that the sum of firing rate remains aproximatively constant.
%80 deg prior
GaussPriorGain = scale2unit * vmPdfs(1:1:360,225,2.7,[]);
drawVoxTuning(GaussPriorGain.*tuningVox135,GaussPriorGain.*tuningVox315,...
    [693 756 617 324],'linestyle','-','color',[1 0 0])

%20 deg prior
GaussPriorGain = scale2unit * vmPdfs(1:1:360,225,8,[]);
drawVoxTuning(GaussPriorGain.*tuningVox135,GaussPriorGain.*tuningVox315,...
    [693 756 617 324],'linestyle','-','color',[1 0.5 0])

%10 deg prior
GaussPriorGain = scale2unit * vmPdfs(1:1:360,225,33,[]);
drawVoxTuning(GaussPriorGain.*tuningVox135,GaussPriorGain.*tuningVox315,...
     [693 756 617 324],'linestyle','-','color',[.3 0.7 .3])

legend('dir 135 - 80 deg prior','dir 315 - 80 deg prior',...
    'dir 135 - 40 deg prior','dir 315 - 40 deg prior',...
    'dir 135 - 20 deg prior','dir 315 - 20 deg prior',...
    'dir 135 - 10 deg prior','dir 315 - 10 deg prior')

%prior mode
hold all
plot([225 225],[0 14],'k--','linesmoothing','on')
text(226,13.5,'Prior mode')

%convention
drawPublishAxis


%-----------
%Competition
%-----------
%Competition is modeled via replacing Bayesian integration by inhibitory
%connections between all readout neurons
figure(2)

%sort trials subjects report the motion direction
%------------------------------------------------
subplot(121)

%voxels tuning
scale2unit = 10;
tuningVox135 = scale2unit*vmPdfs(1:1:360,135,2.4,[]);
tuningVox315 = scale2unit*vmPdfs(1:1:360,315,2.4,[]);

%plot tuning
drawVoxTuning(tuningVox135,tuningVox315,[692 358 617 324],'linestyle','-','color',[0.5 .1 0])
title({'Bayesian integration: stronger priors', 'induce larger shifts in voxels seletivities matching',...
    ' perceptual biases'})

%plot qualitative simulation of voxel selectivity competition readout area 
%when subjects choose the motion direction
GaussPriorGain = scale2unit * vmPdfs(1:1:360,225,2.7,[]);
drawVoxTuning(GaussPriorGain.*tuningVox135,GaussPriorGain.*tuningVox315,...
    [692 358 617 324],'linestyle','-','color',[1 0 0])

%20 deg prior
GaussPriorGain = scale2unit * vmPdfs(1:1:360,225,8,[]);
drawVoxTuning(GaussPriorGain.*tuningVox135,GaussPriorGain.*tuningVox315,...
    [692 358 617 324],'linestyle','-','color',[1 0.5 0])

%10 deg prior
GaussPriorGain = scale2unit * vmPdfs(1:1:360,225,33,[]);
drawVoxTuning(GaussPriorGain.*tuningVox135,GaussPriorGain.*tuningVox315,...
     [692 358 617 324],'linestyle','-','color',[.3 0.7 .3])

legend('dir 135 - 80 deg prior','dir 315 - 80 deg prior',...
    'dir 135 - 40 deg prior','dir 315 - 40 deg prior',...
    'dir 135 - 20 deg prior','dir 315 - 20 deg prior',...
    'dir 135 - 10 deg prior','dir 315 - 10 deg prior')

%prior mode
hold all
plot([225 225],[0 14],'k--','linesmoothing','on')
text(226,13.5,'Prior mode')


%sort trials subjects report the prior mean
%------------------------------------------
subplot(122)

%same but when subjects choose the prior mode




%convention
drawPublishAxis



function drawVoxTuning(tuningVox135,tuningVox315,position,varargin)

set(gcf,'color','w','position',position)
hold all
plot(tuningVox135,'linesmoothing','on',varargin{:})
plot(tuningVox315,'linesmoothing','on',varargin{:})
box off
xlabel('Motion directions (deg)')
ylabel('BOLD (% change)')



