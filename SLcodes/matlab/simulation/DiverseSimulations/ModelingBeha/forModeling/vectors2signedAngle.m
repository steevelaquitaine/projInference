%Calculate the angle formed by two vectors (distanceR2
function angle=vectors2signedAngle(v1,v2)
%Inputs are two vectors' coordinates
%xV1=v1(1)
%yV1=v1(2)
%xV2=v2(1)
%yV2=v2(2)

%e.g.,v1.x=0;v1.y=1;v2.x=1;v2.y=0;
%angle=- (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x,v1.x*v2.x+v1.y*v2.y)
%gives 90 degrees.

%v1.x=1;v1.y=0;v2.x=0;v2.y=1;
%angle=- (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x,v1.x*v2.x+v1.y*v2.y)
%gives - 90 degrees.

%data
xV1=v1(:,1);
yV1=v1(:,2);
xV2=v2(:,1);
yV2=v2(:,2);

%Calculate the angle in degree separating the two vectors
angle=-(180/pi).*atan2(xV1.*yV2-yV1.*xV2,xV1.*xV2+yV1.*yV2);