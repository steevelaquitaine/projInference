
%SLstimsize2visangle.m
%
%  author: Steeve Laquitaine, inspired by 25/9/2009 J Carlin's code
%
% purpose: Convert stimulus size in pixels to degrees visual angle. 
%
%   usage: [StimWidthInVisAngle,StimHeightInVisAngle] =  SLstimsize2visangle(3,3,505,285,385,[1024 768])
%
%       StimWidthInPix  (in deg)
%       StimHeightInPix (in deg)
%       dist2screen     (in mm)
%       screenwidthmm   (in mm)
%       screenHeightmm  (in mm)
%       screenres       (in pixels)
%
%reference: http://imaging.mrc-cbu.cam.ac.uk/imaging/TransformingVisualAngleAndPixelSize


function [StimWidthInVisAngle,StimHeightInVisAngle] = SLstimsize2visangle(StimWidthInPix,StimHeightInPix,dist2screen,screenwidthmm,screenHeightmm,screenres)

%visual angle in deg
visang_radWidth = 2 * atan(screenwidthmm/2/dist2screen);
visang_radHeight= 2 * atan(screenHeightmm/2/dist2screen);

%visual angle in deg
visang_degWidth = visang_radWidth * (180/pi);
visang_degHeight = visang_radHeight * (180/pi);

%visual angle per pix
visang_perpixWidth = visang_degWidth / screenres(1);
visang_perpixHeight = visang_degHeight / screenres(2);

%Stimulus width and Height in visual angle
StimWidthInVisAngle = StimWidthInPix * visang_perpixWidth;
StimHeightInVisAngle = StimHeightInPix * visang_perpixHeight;

