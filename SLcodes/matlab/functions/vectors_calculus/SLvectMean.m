%vector average
%Convert all angles to corresponding points on the unit circle, 
%e.g.,  to (cos(theta),sin(theta)). That is, convert polar coordinates to 
%Cartesian coordinates. Then compute the arithmetic mean of these points. 
%The resulting point will not lie on the unit disk. Convert that point 
%back to polar coordinates.
%v is a directional vector (cartesian coordinates)
%r is the unit circle radius.
%reference: http://en.wikipedia.org/wiki/Mean_of_circular_quantities

function Av=vectMean(v,r)

%cartesian coordinates of the mean
pVectors=sum(v.*r(ones(2,1),:)');

%radian
[Av.rad,Av.deg]=cart2polar(pVectors);


