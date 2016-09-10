      
      %author: steeve laquitaine
        %date: 140426
       %usage: 
               %llhBOLDmn=rand(360,67)
               %[llhBOLDamptSeries,dirTrainedtSeries]=makeBOLDtSeriesByBlockDesign(llhBOLDmn,12)
     %purpose: makes a matrix of BOLD amplitudes (n trials, m voxels)
     %obtained for example from deconvolution into a Block design time 
     %Series. The BOLD response to stimulus presentation is repeated for 
     %the duration of the stimulus.

function [llhBOLDamptSeries,dirTrainedtSeries]=makeBOLDtSeriesByBlockDesign(llhBOLDmn,StimDuration,dirTrained)

%number of voxels
numVoxllh=size(llhBOLDmn,1);

%number of trials
numTrials=numel(dirTrained);

%make a time series of BOLD amplitudes by repeating the amplitude of the 
%BOLD response at each time steps of stimulus presentation. Do that for 
%each voxel.
llhBOLDamptSeries=nan(numTrials*StimDuration,numVoxllh);
for j=1:numVoxllh
    
    %repeat BOLD 
    llhBOLDamptStmp=repmat(llhBOLDmn(j,:),StimDuration,1);
    llhBOLDamptSeries(:,j)=llhBOLDamptStmp(:);
end

%keep track of stimulus at each time steps
dirTrainedtStmp=repmat(dirTrained(1,:),StimDuration,1);
dirTrainedtSeries=dirTrainedtStmp(:);

%print
fprintf('%s \n','(makeBOLDtSeriesByBlockDesign) Time Series of BOLD amplitude...done')

