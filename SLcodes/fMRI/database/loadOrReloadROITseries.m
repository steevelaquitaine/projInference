

%loadOrReloadROITseries.m
%
%
% author: steeve laquitaine
%purpose: load or reload ROITseries fast by saving them in 
%         a directory after each load
%
%usage : 
%             o.myROIname = 'V1';
%             o.sessPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/'		
%             o.myROIpath = '~/data/datafMRI/mlrAnatDB/s0025/mlrROIs/';
%             o.myGroup = 'Averages';
%             o.stckPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/'		
% 
%             %fast if roi data already exist
%             mythisroiloadOrReloadROITseries(v,1,1)
% 
%             slower
%             mythisroi = loadOrReloadROITseries(v,1,1)


function mythisroi = loadOrReloadROITseries(v,thisroi,thissession,o,varargin)

%get args
o.myvarg = varargin;

%re-loadROITseries (very long, not recommended except if there is no
%already existing roi tseries)
if ~slIsInput(o.myvarg,'loadSavedROI')        
    mythisroi = slsaveROItseries(v,o,thisroi,thissession);            
elseif slIsInput(o.myvarg,'loadSavedROI')        
    %or load saved roi tseries (recommended, much faster)    
    [~,mythisroi] = slimportROItseries(o,thissession,thisroi);        
end

%save a session ROI tseries in structured directory
function mythisroi = slsaveROItseries(v,o,thisroi,thissession)

%check if cell or mat
if iscell(v)
    thisv = v{thissession};
else
    thisv = v;
end
if iscell(o.myROIname)
    thisroi = o.myROIname{thisroi};
else
    thisroi = o.myROIname;
end

%load roi tseries
mythisroi = loadROITSeriesFromPath(thisv,thisroi,[],[],['roidir=' o.myROIpath]); %load roi tseries (very slow)

%check cell or mat
if iscell(o.sessPath)
    thissesspath = o.sessPath{thissession};
else
    thissesspath  = o.sessPath;
end

%check roi
if isempty(mythisroi)    
    fprintf('\n')    
    slPrintfStr('slsaveROItseries',' Could not find ROI. Path incorrect ?.....................')    
    fprintf('\n')    
    mrQuit; keyboard    
end

%session path and name
[sPath,fn] = fileparts(thissesspath);

%check dir exists for this session
savedDir = [o.stckPath 'slStckAnalyses/' o.myGroup '/' fn '/roiTseries/' thisroi];
if ~isdir(savedDir)    
    mkdir(savedDir)    
end

%save ROI tseries data
cd(savedDir)
save(['tSeries' o.myGroup fn thisroi date],'mythisroi')
cd ..
