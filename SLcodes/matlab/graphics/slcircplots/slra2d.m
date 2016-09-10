
%slra2d.m
%
%      author: steeve laquitaine
%        date:
%     purpose: convert from radians to degrees
%
%       usage: 
%
%               degrees = slra2d(theta)    %theta in radian
%               
% description: output angle ranges between 1 : 360. It cannot be 0;

function degrees = slra2d(theta)

%When input radians are between 0:2*pi
degrees=(theta/(2*pi))*360;

%if larger than 360 degrees then subtract
%360 degrees
while (sum(degrees>360))
    degrees = degrees - (degrees>360)*360;
end

%if less than 360 degrees then add
%360 degrees
while (sum(degrees<-360))
    degrees = degrees + (degrees<-360)*360;
end

%When radians are signed between -pi:pi.
degrees(degrees<0)=degrees(degrees<0)+360;

%make sure angle ranges between 1:360
degrees(degrees==0)=360;