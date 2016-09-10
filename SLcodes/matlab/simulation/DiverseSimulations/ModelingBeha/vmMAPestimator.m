
%Steeve Laquitaine 130902

%Usage: 
    %di=1:1:360;
    %p is a von Mises pdf.
    %Av=vmMAPestimator(di,p)

%Notes    
    %To be sure that this gives the true mean, I checked with von Mises
    %distributions with known means.
    %di=[1:1:360]';
    %p=vmPdf(di,u,k);

    %note
    %to test a non von Mises distributions (e.g., a uniform pdf) be sure
    %that the range is symetric about the mean. [1:1:359] is good because it
    %is. [1:1:360] is not good. Try to visualize each individual directional
    %vector on the screen to get the idea. If we start with a direction shifted
    %1 degree from 0 whe have to end with a direction shifted -1 degree from 0
    %to ensure perfect symetry around the mean (180);
    %di=[1:1:359]'; p=ones(numel(di),1)/numel(di); Av=vmMAPestimator(di,p)

    %note
    %180 degrees' y cartesian coordinates should be exactly 0, but it isn't
    %exactly at e-10-sthg...

function Av=vmMAPestimator(di,p)
%for test with a von Mises.
%u=180;
%k=42;
%p=vmPdf(di,u,k);


%MAP estimate
%convert directions into vectors
v=polar2cartesian(di,p);

%vectors are averaged.
Av=vectorAverage(v);


%Nested
%vector average
function Av=vectorAverage(v)
%v is a directional vector (cartesian coordinates)
%r is the unit circle radius.

%cartesian of the mean
Av.cart=nanmean(v,1);
x=Av.cart(:,1);
y=Av.cart(:,2);

%radian (arctan(opposite/adjacent))
%Av.rad=atan(y./x);
[Av.rad, Av.deg]=cart2polar(Av.cart);

%polar to cartesian
function v=polar2cartesian(di,r)
%di:angs in degree
%r:radius of the unit circle

%degrees & radians
di2.deg=di;
di2.rad=(di2.deg*pi)/180;
% di2.rad=(di2.deg*3.14)/180;


%cartesian vectors
x=r.*cos(di2.rad);
y=r.*sin(di2.rad);
v=[x y];

%polar to cartesian
function [rad,deg]=cart2polar(v)
%v is the cartesian vector

%radian (arctan(opposite/adjacent))
x=v(:,1);
y=v(:,2);
rad=atan(y./x);

%degrees
deg=ra2d(rad);

%adjust each angle according to his quadrant (in degree) otherwise atan
%gives wrong values for negative x and y.
for i=1:numel(x)
    if x(i)>=0 && y(i)>=0              %(quadrant 1)
        deg(i)=deg(i);
    elseif x(i)<0                      %(quadrant 2 & 3)
        deg(i)=deg(i)+180;
    elseif x(i)>=0 && y(i)<0           %(quadrant 4)
        deg(i)=deg(i)+360;
    end
end

%degrees to radians
function radians=de2r(ang)
radians=(ang/360)*2*pi;

%convert radians to degrees
function degrees=ra2d(theta)
degrees=(theta/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end

%von Mises pdf
function mPdfs=vmPdf(x,u,k)
%x: space
%u: mean
%k: precision or concentration parameter
w=1;

%degrees to radian
u=de2r(u);
x=de2r(x);

%pdf(sum to 1)
mPdfs=exp(k.*cos(w*(x-u)))/2*pi.*besseli(0,k);
mPdfs=mPdfs/sum(mPdfs);




