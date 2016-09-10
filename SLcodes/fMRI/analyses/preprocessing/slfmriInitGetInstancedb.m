

%slfmriInitGetInstancedb.m
%
% author: steeve laquitaine
%purpose: initialize parameters to run "slfmriGetInstancedb"
%usage:
%
%       [o,sessPath] = slfmriInitGetInstancedb('V1','ClassifStckSessmyRandomCoh_V1_23-Feb-2016.mat');
%

function [o,sessPath] = slfmriInitGetInstancedb(roiname,filename)

o.roiname = roiname;
% o.dataPath = ['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/' o.roiname '/'];
% sessPath{1} = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
% sessPath{2} = '~/data/datafMRI/sltaskdotdirfmri05/s002520150923';
% sessPath{3} = '~/data/datafMRI/sltaskdotdirfmri05/s002520150925';
% o.filename = filename;

if isfield(o,'sessPath')
    o = rmfield(o,'sessPath');
end


