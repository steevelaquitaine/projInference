% mglMakeGaussian.m
%
%      usage: mglMakeGaussian(width,height,sdx,sdy,<xCenter>,<yCenter>,<xDeg2pix>,<yDeg2pix>)
%         by: justin gardner
%       date: 09/14/06
%    purpose: make a 2D gaussian, useful for making gabors for
%             instance. Set mglVisualAngleCoordinates before
%             using.
%
%             width, height, sdx and sdy are in degrees of visual angle
%
%             xcenter and ycenter are optional arguments in degrees of visual angle
%             and default to 0,0
%
%             xDeg2pix and yDeg2pix are optional arguments that specify the
%             number of pixels per visual angle in the x and y dimension, respectively.
%             If not specified, these values are derived from the open mgl screen (make
%             sure you set mglVisualAngleCoordinates).
%       e.g.:
%
%        e.g., gaussianRing = mglMakeGaussianRingRentzi(5,5,0.5,1/64,0.25,'expDecay')
%
%
%inputs:
%       type: 'expDecay' or 'linearDecay'
%
% mglOpen;
% mglVisualAngleCoordinates(57,[16 12]);
% mglClearScreen(0.5);
% grating = mglMakeGrating(10,10,1.5,45,0);
% gaussian = mglMakeGaussian(10,10,1,1); 
% gabor = 255*(grating.*gaussian+1)/2;
% tex = mglCreateTexture(gabor);
% mglBltTexture(tex,[0 0]);
% mglFlush;



function gaussianRing = mglMakeGaussianRingRentzi(width,height,innerCircle,outerCircle,stdRatio,type)
%stdRatio is the std with respect to the number of elements contained in
%the lineGaussian. for example stdRatio = 1/4 indicates std = LengthGaussian/4
% innerCircle = 1/2;
% outerCircle = 1/64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHECK ARGUMENTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xCenter = 0; yCenter = 0;

if isempty(mglGetParam('xDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
end
xDeg2pix = mglGetParam('xDeviceToPixels');


if isempty(mglGetParam('yDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
end
yDeg2pix = mglGetParam('yDeviceToPixels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get size in pixels
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
%if the width or height in Pixels are not an odd number we make them an odd number
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

% get a grid of x and y coordinates that has 
% the correct number of pixels
x = -width/2:width/(widthPixels-1):width/2;
y = -height/2:height/(heightPixels-1):height/2;
[xMesh,yMesh] = meshgrid(x,y);

% compute gaussian window
m = [];
sdx = width/2;
sdy = height/2;
m = exp(-(((xMesh-xCenter).^2)/(2*(sdx^2))+((yMesh-yCenter).^2)/(2*(sdy^2))));
% clamp small values to 0 so that we fade completely to gray.
m(m(:)<0.01) = 0;


%%%%%%%%%%%%%%%%%%%%%%  MAKE INNER AND OUTER MASK  %%%%%%%%%%%%%%%%%%%%%%%%%
partOutMask= 1 - 1*(m>=exp(-outerCircle)); %The part Out of the big circle is 0, the part in is 1
partInMask = 1 - 1*(m<=exp(-innerCircle)); %The part in the small circle is 0, the part out is 1
fullMask = partOutMask + partInMask -1; %Full mask, where you keep it is 1, when not it is 0
gaussianMasked = fullMask.*m;
centerX = ceil(size(gaussianMasked,2)/2);
centerY = ceil(size(gaussianMasked,1)/2);

%BELOW NESTED IN THE IF STRUCTURE, WE FIND THE INDSTART AND INDENDVVARIABLES 
%THAT WILL BE USED TO FIND THE SIZE OF THE GAUSSIAN WE AREVGOING TO USE
%width is always wider, but just in case
if centerX > centerY   %if width is greater than height
    
    indStart = centerX;
    flag = 0;
    while flag == 0
        % find the starting nonzero poing on the horizontal axis of the figure
        if gaussianMasked(centerY,indStart) == 0
            indStart = indStart +1;
        else
            flag = 1;
        end
    end
    
    indEnd = size(gaussianMasked,2);
    flag = 0;
    while flag == 0
        % find the starting nonzero poing on the horizontal axis of the figure
        if gaussianMasked(centerY,indEnd) == 0
            indEnd = indEnd - 1;
        else
            flag = 1;
        end
    end
    
else
    sprintf('problem with the dimensions. Height is greater than width')
    return
end

%set the matrix with the indices that we are going to use for the Gaussian
valuesToCompare = gaussianMasked(centerY,indStart:indEnd);
indMatrix = zeros(size(gaussianMasked,1),size(gaussianMasked,2));
for i = 1: size(gaussianMasked,1) % height
    for j = 1:size(gaussianMasked,2) %width
        
        if gaussianMasked(i,j)~= 0
            [val, indMatrix(i,j)] = min(abs(valuesToCompare - gaussianMasked(i,j)));
        end
    end
end

keyboard

%%MAKE THE LINE GAUSSIAN WE ARE GOING TO USE AROUND THE RING
lengthGaussian = indEnd - indStart +1;
xx = 1:1:lengthGaussian;
lineCenter = round(lengthGaussian/2);
lineStd = lengthGaussian*stdRatio;

%%%%%%% case exponential decay
if strcmp(type,'expDecay')    
    gaussianLine = exp(-(((xx-lineCenter).^2)/(2*(lineStd^2))));
end

%%%%%%% case linear decay
if strcmp(type,'linearDecay')   
    %the larger the stdRatio the larger the area of linear decay
    
    %margin area where linear decay is not aplied
    margin = max([round((1-stdRatio)*lengthGaussian) 1]);
    if margin>lengthGaussian/2 +1
        margin = round(lengthGaussian/2) - 1;
    end
    %area where actual linear decay is applied
    range = margin:lengthGaussian-margin;
    gaussianLine = abs((xx(range)-lineCenter));
    %scale between 0 and 1
    %gaussianLine = 1 - (gaussianLine-min(gaussianLine(:)))/(max(gaussianLine(:))-min(gaussianLine(:)));
    gaussianLine = (gaussianLine-min(gaussianLine(:)))/(max(gaussianLine(:))-min(gaussianLine(:)));
    %pad margin with zeros
    gaussianLine1 = [zeros(1,margin) gaussianLine];
    gaussianLine  = [gaussianLine1 zeros(1,lengthGaussian-length(gaussianLine1))];  
end


%Make the GaussianRing
gaussianRing = zeros(size(gaussianMasked,1),size(gaussianMasked,2));
for i = 1:lengthGaussian
    
    indG = find(indMatrix == i);
    gaussianRing(indG) = gaussianLine(i);
    
end

% imshow(gaussianRing)
