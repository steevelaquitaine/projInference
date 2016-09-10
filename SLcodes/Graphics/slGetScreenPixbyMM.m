
%slGetScreenPixbyMM.m
%
%author: steeve laquitaine
%  date: 151203
%purpose: get pixels by mm on screen
%
%usage:
%
%       PixBymm = slGetScreenPixbyMM
%
%require mgl (justin gardner) library


function PixBymm = slGetScreenPixbyMM(screenNumber)


if screenNumber==0
    screenNumber=1;
end

%convert distance (in pix) to visual angle (deg)
dinfo = mglDescribeDisplays;
%get called screen
dinfo = dinfo(screenNumber);
%screen true width in mm
widthmm = dinfo.screenSizeMM(1);
%screen true height in mm
heightmm = dinfo.screenSizeMM(2);
%screen width in pix
widthpix = dinfo.screenSizePixel(1);
%screen height in pix
heightpix = dinfo.screenSizePixel(2);
%mm per pix
PixBymm = min([widthpix/widthmm heightpix/heightmm]);
