

%slmakeDisk.m
%
%
% author: steeve laquitaine
%purpose: draw a disk



function m = slmakeDisk(radius,xCenter,yCenter)

%check arguments
m = [];

if ieNotDefined('xCenter'),xCenter = 0;end
if ieNotDefined('yCenter'),yCenter = 0;end

%defaults for xDeg2pix
if ieNotDefined('xDeg2pix')
  if isempty(mglGetParam('xDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
  end
  xDeg2pix = mglGetParam('xDeviceToPixels');
end

%defaults for yDeg2pix
if ieNotDefined('yDeg2pix')
  if isempty(mglGetParam('yDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
  end
  yDeg2pix = mglGetParam('yDeviceToPixels');
end

%get size in pixels
widthPixels = round(radius*xDeg2pix);
heightPixels = round(radius*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

%get a grid of x and y coordinates that has 
%the correct # of pixels
x = -radius : radius/(widthPixels-1) : radius;
y = -radius : radius/(heightPixels-1) : radius;
[xMesh,yMesh] = meshgrid(x,y);

%compute disk
m = sqrt((xMesh-xCenter).^2 + (yMesh-yCenter).^2) < radius;
