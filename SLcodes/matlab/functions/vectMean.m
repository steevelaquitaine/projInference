 %Author: Steeve Laquitaine
   %date: last update 140107
  %usage: [deg,rad,pVectors]=vectMean([90;0], [1 1])
         %or [deg,rad,pVectors]=vectMean([90;0], [1 1])
         %v: directional vector (cartesian coordinates)
         %r: unit circle radius.
%purpose: calculate vector average. Convert all angles to corresponding 
%points on the unit circle, e.g.,  to (cos(theta),sin(theta)). That is, 
%convert polar coordinates to cartesian coordinates. Then compute the 
%arithmetic mean of these points. The resulting point will not lie on 
%the unit disk. Convert that point back to polar coordinates.

   %note: 
   %.vector average of 270 and 90 produces 180. Even if counterintuitive
   %it is correct.
   
   %.Vector mean depends on the angles' sequence order. I think this should
   %not be the case, but it's because of round off errors for small
   %numbers. For now I don't know how to solve this...
   %deg1=vectMean([0;90;180;270],[1 1 1 1])
   %deg2=vectMean([270;0;90;180],[1 1 1 1])
   
%reference: http://en.wikipedia.org/wiki/Mean_of_circular_quantities
%reference: http://stackoverflow.com/questions/1813483/averaging-angles-again

function [deg,rad,pVectors]=vectMean(v,r)

%if polar angles are input convert to coordinates
if size(v,2)==1
    v=polar2cartesian(v,1);
end

%mean coordinates
weight=r(ones(2,1),:)';
pVectors=sum(v.*weight,1);

%radian
[rad,deg]=cart2polar(pVectors);