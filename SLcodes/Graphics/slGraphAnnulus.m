% slGraphAnnulus.m
%
%     author: steeve laquitaine
%       date: 16/02/06
%    purpose: draw an annulus
%            
%      usage:
%
%        e.g., use as a standalone
%
%               annulus = slGraphAnnulus(20,20,1.5,2,0.5);
%               figure; imagesc(annulus)
%
%       ?e.g.,?use with mgl libary
%
%               mglOpen(0);
%               mglVisualAngleCoordinates(57,[16 12]);
%               mglClearScreen(0.5);
%               grating = mglMakeGrating(20,20,1.5,45,0);
%               annulus = slGraphAnnulus(20,20,1.5,2,0.5);
%               gabor = 255*(grating.*annulus+1)/2;
%               tex = mglCreateTexture(gabor);
%               mglBltTexture(tex,[0 0]);
%               mglFlush;
%
%
% Inputs:
%               width: Screen width (degrees)
%              height: Screen height (degrees)
%         innerRadius: distance from the screen center to the inner edge 
%                      of the annulus (degrees)      
%         outerRadius: distance from the screen center to the outer edge 
%                      of the annulus (degrees)      
%       decayDistance: annulus linearly decays from inner/outer edge up to 
%                     "decayDistance" inward/outward. i.e.,  at inner/outer egde
%                      radius -/+ "decayDistance" annulus has completely
%                      disappeared (0 values) to leave place to the
%                      background
%                      (degrees)
%
%reference: inspired from code from Ilias Rentzeperis and Justin Gardner



function AnnulusMask = slGraphAnnulus(width,height,innerRadius,outerCircle,decayDistance)

%convert from degrees to pixels
if isempty(mglGetParam('xDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
end
if isempty(mglGetParam('yDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
end
xDeg2pix = mglGetParam('xDeviceToPixels');
yDeg2pix = mglGetParam('yDeviceToPixels');

%get annulus size in pixels
%If they are not, make width and height an odd number
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

%get a grid of x and y coordinates that has 
%the correct number of pixels (grid values in degrees)
x = -width/2  : width/(widthPixels-1)   : width/2;
y = -height/2 : height/(heightPixels-1) : height/2;
[xMesh,yMesh] = meshgrid(x,y);


%%%%%%%%%%%%%%%%%%%%%%  MAKE INNER AND OUTER MASK  %%%%%%%%%%%%%%%%%%%%%%%%%
%screen center
center = [0 0];

%coord to matrix center (degrees)
x2c = xMesh - center(2);       %y to center
y2c = -(yMesh - center(1));    %x to center 

%create distance-to-center map in degrees
r = sqrt(y2c.^2 + x2c.^2); %radius (Pythagoras)

%create annulus mask
AnnulusMask = ones(heightPixels,widthPixels);
AnnulusMask(r<innerRadius)=0;
AnnulusMask(r>outerCircle)=0;

%checking
checkOuter = r>outerCircle;
if ~any(checkOuter==1)
    fprintf('(slGraphAnnulus) Annulus outer edgde is too large for the screen. Reduce it.')
    mglClose; keyboard
end

%implement linear decaying to 0 over the inner .5 degree radius
%inward linear decaying mask
DecayArea = zeros(heightPixels,widthPixels);
DecayAreaIn = r>(innerRadius-decayDistance) & r<innerRadius;
DecayArea(DecayAreaIn) = r(DecayAreaIn);
x = DecayArea(DecayAreaIn);
%from out to in
scaledx = (x - min(x(:)))/(max(x(:)) - min(x(:)));
DecayArea(DecayAreaIn) = scaledx;

%outward linear decaying mask
DecayAreaOut = r<(outerCircle+decayDistance) & r>outerCircle;
DecayArea(DecayAreaOut) = r(DecayAreaOut);
x = DecayArea(DecayAreaOut);
%1 - (...) because from in to out
scaledx = 1 - ((x - min(x(:)))/(max(x(:)) - min(x(:))));
DecayArea(DecayAreaOut) = scaledx;

%Linearly decaying onward and outward Annulus
AnnulusMask = AnnulusMask + DecayArea;
