
%slGraphWedges2.m
%
% author: steeve laquitaine
%   date: 151216
%purpose: create a matrix of bright oriented wedge area (1) on a black
%         background (0)
%
%  usage:
%
%ex1:
%       Wdg = slGraphWedges2(400,400,100,315,20)
%       imagesc(Wdg);colormap('gray')
%
%ex2:
%       Wdg = slGraphWedges2(400,400,100,315,20)
%       mytex = mglCreateTexture(flipud(Wdg));      
%      
%
%
%Description: 
%
%       H: Matrix (image) Height in pix
%       W: Matrix (image) Width in pix
%      wr: wedge radius in pix
%      th: wedge orientation (deg)
%   tspan: wegde angle span (deg)
%
%     Wdg: a matrix of 1 (wedge) and 0 (backg)
%          imshow(Wdg) or imagesc(Wdg); colormap('gray')
%
%
%   Use some trigonometrics to create wedges as overlap of
%   angles and distance maps
%   tan(ang)=y/x; x=r*cos(ang); y=r*sin(ang); 





function Wdg = slGraphWedges2(H,W,wr,th,tspan)

%2D space point coordinates
[x,y] = meshgrid(1:W,1:H);
center = [H/2 W/2];

%coord to matrix center
x2c = x - center(2);       %y to center
y2c = -(y - center(1));    %x to center 

%create distance map (to center)
%in matrix elmts
r = sqrt(y2c.^2 + x2c.^2); %radius (Pythagoras)
 
%create angle map (to center)
ang = atan2(y2c,x2c);      %rad
angd = round(SLra2d(ang)); %deg

%create oriented wedge
%as the map intersection. 
%Wedge spans tspan around th angle
th2AngDcw = SLra2d(SLde2r(th,1)-ang); %map of angles distant by +tspan/2
th2AngDccw = SLra2d(ang-SLde2r(th,1));%by -th/2 map
fllWdg = (th2AngDcw<=(tspan/2) | th2AngDccw<=(tspan/2)); %full wedge

%wedge is intersection of maps
Wdg = zeros(H,W);
Wdg(r<wr & fllWdg==1)=1;
%imshow(Wdg)
%imagesc(Wdg); colormap('gray')

