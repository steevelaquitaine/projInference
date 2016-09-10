
%slsimfMRIdatabase.m
%
% author: steeve laquitaine
%purpose: simulate fMRI database with equiprobable distribution of 5
%         voxel selectivities (500 voxels, 10 trials, 3 variables)
%
%usage:
%
%       d = slsimfMRIdatabase
%
%Inputs
%       no input : 5 voxel selectivities equiprobably distributed
%       'random': database is filled with random values
%       'vmDist': direction are von mises dist. across trials
%
%Descriptions:
%
%   simulates 3 variables: "mySwitch","myRandomDir","myRandomCoh"
%   10 trials (instances) and 500 voxels
%
%
function d = slsimfMRIdatabase(varargin)

%% simulate dataset with equiprobable distribution of direction preferences
%dataset for switch=1, 5 trials, each direction, one coh, 5 select. equip
%distributed
d.mySwitch    = [1 1 1 1 1];
d.myRandomDir = [15 85 155 225 295];
d.myRandomCoh = [.08 .08 .08 .08 .08];
sel = repmat([15 85 155 225 295],1,100);
ks = 4*ones(1,length(sel));
d.instances = vmPdfs(d.myRandomDir,sel,ks,'norm');

%identical dataset for switch 2
d.mySwitch = [d.mySwitch 2 2 2 2 2]';
d.myRandomDir = [d.myRandomDir 15 85 155 225 295]';
d.myRandomCoh = [d.myRandomCoh .08 .08 .08 .08 .08]';
d.instances = [d.instances;d.instances];
d.stimvols = nan(length(d.myRandomDir),1);

%case random values
if strcmp(varargin,'random')
    d.instances = rand(size(d.instances));
end

%case von mises distribution of directions across trials
if strcmp(varargin,'vmDist')
    d.mySwitch    = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    d.myRandomDir = [15 85 85 155 155 155 225 225 225 225 225 225 295 295 295];
    d.myRandomCoh = [.08 .08 .08 .08 .08 .08 .08 .08 .08 .08 .08 .08 .08 .08 .08];
    sel = repmat([15 85 155 225 295],1,100);
    ks = 4*ones(1,length(sel));
    d.instances = vmPdfs(d.myRandomDir,sel,ks,'norm');  
    
    %identical dataset for switch 2
    d.mySwitch = [d.mySwitch d.mySwitch*2]';
    d.myRandomDir = [d.myRandomDir d.myRandomDir]';
    d.myRandomCoh = [d.myRandomCoh d.myRandomCoh]';
    d.instances = [d.instances;d.instances];
    d.stimvols = nan(length(d.myRandomDir),1);    
end




