

%slmakeAnnulus.m
%
%
% author: steeve laquitaine
%   date: 151210
%purpose: draw an annulus


function m = slmakeAnnulus(radiusOut,xCenter,yCenter)

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
widthPixels = round(radiusOut*xDeg2pix);
heightPixels = round(radiusOut*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

%get a grid of x and y coordinates that has 
%the correct # of pixels
x = -radiusOut : radiusOut/(widthPixels-1) : radiusOut;
y = -radiusOut : radiusOut/(heightPixels-1) : radiusOut;
[xMesh,yMesh] = meshgrid(x,y);

%compute annulus
maxRad = sqrt((xMesh-xCenter).^2 + (yMesh-yCenter).^2) < radiusOut;

%inner radius is .5 degrees inward
radiusIn = radiusOut-0.2;
minRad = sqrt((xMesh-xCenter).^2 + (yMesh-yCenter).^2) > radiusIn;
m = zeros(size(maxRad));
m(minRad==1 & maxRad == 1) = 1;
m(minRad==0 & maxRad == 0) = 0;


