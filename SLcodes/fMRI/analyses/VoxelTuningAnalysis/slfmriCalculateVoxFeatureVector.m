
%slfmriCalculateVoxFeatureVector.m
%
%
%author : steeve laquitaine
%purpose: calculate voxels feature vectors by using Georgopoulos et al.'s
%         definition of population vector for a voxel responses to the
%         different states of a stimulus feature (e.g., different motion 
%         directions)
%
%
%y, polar angles, and radius must be column vectors


function [degmean,coordmean] = slfmriCalculateVoxFeatureVector(y,r)

%cartesians
rads = y*pi/180;
xx = round(r.*cos(rads).*10^4)./10^4;
yy = round(r.*sin(rads).*10^4)./10^4;

%cartesian weighted means
coordmean = r'*[xx yy];
degmean   = SLgetangle(coordmean(:,1),coordmean(:,2));