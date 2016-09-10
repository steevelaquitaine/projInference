
 %author: steeve laquitaine
   %date: 131128
  %usage: [rad,deg]=vectMean([0 1; 1 0],[0.75 0.25])
  %you can draw the vectors to with drawVectors(Av.deg,1)
   %purpose: calculate vector averages
         %convert all angles to corresponding points on the unit circle, 
         %e.g., to (cos(theta),sin(theta)). That is, convert polar 
         %coordinates to cartesian coordinates. Then compute the 
         %arithmetic mean of these points. The resulting point will not 
         %lie on the unit disk. Convert that point back to polar 
         %coordinates.
%v is a directional vector (cartesian coordinates)
%r is the unit circle radius.
%reference: http://en.wikipedia.org/wiki/Mean_of_circular_quantities

function [deg,rad]=vectMean(v,r)

%cartesian coordinates of the mean
pVectors=sum(v.*r(ones(2,1),:)');

%radian & degrees
[rad,deg]=cart2polar(pVectors);


