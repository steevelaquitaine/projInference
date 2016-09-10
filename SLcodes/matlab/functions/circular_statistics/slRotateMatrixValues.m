
%slRotateMatrixValues.m
%
%
% author: steeve laquitaine
%   date: 150803
%purpose: rotate radian values of a 
%         matrix x by r clockwisely 
%         (r>0) and ccw (r<0)
%
%
%usage :
%
% 		x = circ_vmrnd(-3.14*ones(5,4),2.7);
% 		rx = slRotateMatrixValues(x,-pi,'dispOn');
%
%		%check rotation
%		for i = 1 : size(x,2) 
%			ro(:,i) = SLvectors2signedAngle(SLra2d(rx(1,:)),SLra2d(x(1,:)),'polar')
%		end
%
%inputs:
%		
%		r : rotation in radian (can be signed)
%		x : matrix of angles in radians
%
%
%ccw rotation matrix: 
%
%			| cos(r)   sin(r) |
%			| -sin(r)   cos(r) |
%
%is multiplied with the corresponding cartesian coordinates of the radian 
%values of each element of x .
%
%note: r = 0 produces no rotation

function rx = slRotateMatrixValues(x,r,dispOpt)

%radian to cart
coo = SLpolar2cartesian(x(:),1,'radian');

%rotate by r
rotMat = [ cos(r) sin(r); 
		  -sin(r) cos(r)];

%Rotation is a matrix multiplication which equivals this for loop
%		for i = 1 : 20; test(i,:) = rotMat*coo(i,:)'; end
%The values of each element of the matrix x are rotated by r
newCoo = (rotMat*coo')';

%back to radians
rxVec = SLcart2polar(newCoo);

%reshape like size(x)
rx = reshape(rxVec,size(x,1),size(x,2));


%plot
if strcmp(dispOpt,'dispOn')
	%plot density of original values directions
    u = 2*pi/360;
	hist(x(:),-pi : u : pi-u)
	hold on; hist(rx(:),-pi : u : pi-u)
end