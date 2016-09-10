

%slfmriGetOverlayData.m
%
%author: steeve laquitaine
%
% usage: 
%
%		myROI = slfmriGetOverlayData()

function myROI = slfmriGetOverlayData()


%load ROI
myROI = loadROITSeries(getMLRView);

%converting the scan coordinates of the ROI to linear coordinates
scanDims = viewGet(getMLRView,'scanDims');

myROI.linearScanCoords = sub2ind(scanDims,myROI.scanCoords(1,:),myROI.scanCoords(2,:),myROI.scanCoords(3,:));

%get r2 map
r2 = viewGet(getMLRView,'overlayData');

%select ROI's r2s
myROI.r2 = r2(myROI.linearScanCoords);


%------------
%plot density
%------------
hist(myROI.r2)
xlabel('R-squared')
ylabel('Count')
box off
set(gcf,'color','w')

