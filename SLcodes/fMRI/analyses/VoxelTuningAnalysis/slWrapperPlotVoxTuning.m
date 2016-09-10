
%slWrapperPlotVoxTuning.m

%% initialize analysis
slfmriInitAnalysisTaskDotDirfMRI05


for voxnum = 26%:299
    st = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch','myRandomCoh=0.12',d,'fminsearch');
    vline(225,'k:')
    pause(0.3)
end