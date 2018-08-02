
%SLvectors2signedAngle.m
%
%    author: steeve laquitaine
%      date: 140529
%   purpose: calculate the angle distance between two vectors in degrees
%            when the 2nd angle is clockwise to the 1st the angle is positive
%
%     usage:
%           angle = SLvectors2signedAngle(90,45,'polar')
%           angle = SLvectors2signedAngle([0 1],[1 0],'cartesian')
%
% varargin:
%
%           'polar'
%           'radian'
%           'cartesian'
%
%
%  Checking:
%       e.g.,v1.x=0;v1.y=1;v2.x=1;v2.y=0;
%       angle=- (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x,v1.x*v2.x+v1.y*v2.y)
%       gives 90 degrees.
%
%       v1.x=1;v1.y=0;v2.x=0;v2.y=1;
%       angle=- (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x,v1.x*v2.x+v1.y*v2.y)
%       gives - 90 degrees.
%
% v1 and v2 should be column vectors if they are in degrees

function angle = SLvectors2signedAngle(v1,v2,varargin)

%check varargin
if isempty(varargin)
    fprintf('SLvectors2signedAngle) Please input angle unit \n')
    fprintf('SLvectors2signedAngle) Your choices are: \n')
    fprintf('SLvectors2signedAngle) - "polar" \n')
    fprintf('SLvectors2signedAngle) - "radian" \n')
    fprintf('SLvectors2signedAngle) - "cartesian" \n')
    keyboard
end

varg = varargin{:};

%make cartesian
if sum(strcmp(varg,'polar')) || sum(strcmp(varg,'radian'))
    
    v1 = SLmakeColumn(v1);
    v2 = SLmakeColumn(v2);
    v1 = SLpolar2cartesian(v1,1,varg); %make cartesian
    v2 = SLpolar2cartesian(v2,1,varg);
    
elseif sum(strcmp(varg,'cartesian'))
    
end

%get cartesian coordinates
xV1 = v1(:,1);
yV1 = v1(:,2);
xV2 = v2(:,1);
yV2 = v2(:,2);

%Calculate the angle in degrees separating the two vectors
angle = -(180/pi).*atan2(xV1.*yV2-yV1.*xV2,xV1.*xV2+yV1.*yV2);