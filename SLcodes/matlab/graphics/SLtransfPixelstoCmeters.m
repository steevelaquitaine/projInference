
%SLtransfPixelstoCmeters.m
%
%   author: steeve laquitaine
%    usage: 
%
%       sizecm = SLtransfPixelstoCmeters
%
%reference: http://stackoverflow.com/questions/3007773/convert-pixel-to-cm
%
%Need to be checked not sure works properly...


function sizecm = SLtransfPixelstoCmeters
          
%find resolution of your display (cm)
pixperinch=get(0,'ScreenPixelsPerInch'); 
inchperpixels=1/pixperinch;
cmperpix=inchperpixels*2.54;

%how many pixels
sizePix=get(0,'ScreenSize'); 

%size cm
sizecm=sizePix([3,4])*cmperpix;            
