
%SLcircCcwVect2Angle.m
%
%    author: steeve laquitaine
%      date: 140807
%   purpose: calculate the angle distance between two vectors in degrees
%            when the second angle is more clockwise the angle is positive
%     usage: 
%           ccwAngle = SLcircCcwVect2Angle(360,270)

function ccwAngle = SLcircCcwVect2Angle(v1,v2)

ccwAngle = vectors2signedAngle(v2,v1);
ccwAngle(ccwAngle<0) = ccwAngle(ccwAngle<0)+360;