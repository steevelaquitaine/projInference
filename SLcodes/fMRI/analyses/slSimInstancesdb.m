
%slSimInstancesdb.m
%
%
% author: steeve laquitaine
%   date: 160121
%purpose: simulate N instances by M dimensions (e.g., voxels) instances
%         for a population of 360 voxels (cols). Voxels selectivities are 
%         circular values and range from 1:1:360 degrees.
%         Each row instance is the response of the 360 voxels to the variable
%         value displayed at this row.
%
%usage: 
%       instances = slSimInstancesdb(d.myRandomDir,1:1:360,10)
%
%inputs
%     var0TSeries : time series of variable 0 circular values (e.g., directions)
%                   e.g., [0, 1, 2 ....]
%       voxSelect : voxels selectivities for variable 0 (e.g., directions, orientation)
%            gain : scaling amplitude of voxel response

function instances = slSimInstancesdb(var0TSeries,voxSelect,gain)

%create 360 example var0 (e.g., directions) tuning curves
%tunings(1,:) is the response to var0(1) of voxels 1 to 360 (with pref 1:1:360);
var0space = 1:1:360;
tunings = gain * vmPdfs(var0space,voxSelect,4,'norm');

%get voxels response for each var0 in the timeseries
for i = 1 : length(var0TSeries)
    
    %this value of var0
    thisVar0 = var0TSeries(i);
        
    %instances
    instances(i,:) = tunings(thisVar0,:);
end

