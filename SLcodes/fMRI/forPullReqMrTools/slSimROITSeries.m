
%slSimROITSeries.m
%
%    Author: steeve laquitaine
%      date: 151008
%   purpose: simulate ROI Tseries
%            The BOLD signal is simulated for each trial at event onset 
%            with a different amplitude for each condition.
%            e.g., all condition 1 trials (e.g., trial 1,10,20) will have a 
%            BOLD value of 1; condition 2 trials a value of 2, etc ....
%            roi.tSeries is a Nvoxels by Nvols matrix. All voxels have the
%            same tSeries.
%           
%            This might be useful if you want to check the validity of
%            classification or event-related analyses or any analysis that 
%            relies on roi tSeries.
%
%     usage:
%
%           v = mrLoadRet([]);    
%
%           %get stimvol for my variable
%           o.stimvol = getStimvol(v,'myRandomDir','taskNum=2','phaseNum=1','segmentNum=2');
%           
%           %load ROI tSeries
%           roi = loadROITSeries(v,'V1'); 
%            
%           %replace roi tSeries by simulations
%           roi = slSimROITSeries(roi,o);
%
%           %plot voxel 1 tSeries with red tagged condition 1 trials
%           plot(roi.tSeries(1,:))        
%           vline(o.stimvol{1},'--r') 
%           xlabel('vols'); ylabel('Bold amplitude'); legend('tSeries','Condition 1 trials')


%simulate ROI Tseries to check the code
function roi = slSimROITSeries(roi,o)

%define a period after event onset
[d2,pos] = sort([o.stimvol{:}]);
d2st  = d2(1:end-1);
d2end = d2(2:end);

%simulate TSeries values for each condition
%value C1 = 1
%value C2 = 2 etc...
roi.tSeries = zeros(size(roi.tSeries));
TSeriesBold = [];
gain = 2;
for i = 1: length(o.stimvol)
    TSeriesBold = [TSeriesBold gain*i*ones(1,length(o.stimvol{i}))];
end
TSeriesBoldsort = TSeriesBold(pos);

%update tSeries with condition simulated Bold values
%all voxels have the same tSeries
for i = 1 : length(d2)-1
    roi.tSeries(:,d2st(i) : d2end(i)-1) = TSeriesBoldsort(i);
end

%fill up 1st trials vols before 1st event
roi.tSeries(:,1:d2st(1)-1) = TSeriesBoldsort(1);

%update last trial value
roi.tSeries(:,d2end(i):d2end(i)+length(d2st(i):d2end(i)-1)) = TSeriesBoldsort(i+1);

%add small noise
%roi.tSeries = roi.tSeries + normrnd(0,0.1,size(roi.tSeries));

fprintf('%s \n','(slSimROITSeries) Roi tSeries have been simulated.')