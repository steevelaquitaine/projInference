 %author: steeve laquitaine
   %date: 140208
%purpose: convert visual angles to pixels
  %usage: size_in_deg=convertVisualAngles2pixels(28.5,50.5,1280,0.07)
            %where:
            %h=28.5 %Monitor height in cm
            %d=50.5 %Distance between monitor and participant in cm
            %r=1280 %Vertical resolution of the monitor
            %size_in_deg=3; %The stimulus size in degrees

%reference:
%http://osdoc.cogsci.nl/miscellaneous/visual-angle/#convert-pixels-to-visual-degrees

function size_in_px=convertVisualAngles2pixels(h,d,r,size_in_deg)

%Calculate the number of degrees that correspond to a single pixel. This will
%generally be a very small value, something like 0.03.
deg_per_px=ra2d(atan2(.5*h, d))/(.5*r);

%Calculate the size of the stimulus in degrees
size_in_px=size_in_deg/deg_per_px;