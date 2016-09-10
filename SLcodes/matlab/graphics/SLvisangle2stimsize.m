
%SLvisangle2stimsize.m
%
%  author: Steeve Laquitaine, inspired by 25/9/2009 J Carlin's code
%
% purpose: Provides x,y size in pixels to produce a given size in visual angle. 
%
%   usage: [StimWidthInPix,StimHeightInPix] =  SLvisangle2stimsize(2.5,2.5,505,300,400,[1920 1200])
%
%       StimWidthInVisangle (in deg)
%       StimHeightInVisangle (in deg)
%       dist2screen   (in mm)
%       screenwidthmm   (in mm)
%       screenheightmm   (in mm)
%       screenres     (in pixels)
%
%reference: http://imaging.mrc-cbu.cam.ac.uk/imaging/TransformingVisualAngleAndPixelSize


function [StimWidthInPix,StimHeightInPix] = SLvisangle2stimsize(StimWidthInVisangle,StimHeightInVisangle,dist2screenmm,screenwidthmm,screenHeightmm,screenres)

%visual angle in deg
visang_radWidth = 2 * atan(screenwidthmm/2/dist2screenmm);
visang_radHeight= 2 * atan(screenHeightmm/2/dist2screenmm);

%visual angle in deg
visang_degWidth = visang_radWidth * (180/pi);
visang_degHeight = visang_radHeight * (180/pi);

%pixels per visual angle
pix_pervisangWidth = screenres(1) / visang_degWidth;
pix_pervisangHeight = screenres(2) / visang_degHeight;

%Stimulus width and Height in pixels
StimWidthInPix = round(StimWidthInVisangle * pix_pervisangWidth);
StimHeightInPix= round(StimHeightInVisangle * pix_pervisangHeight);
