
%slimportROItseries
%
%
% author: steeve laquitaine
%purpose: Import saved ROI tSeries for speed


%load already saved ROI tseries in workspace
function [thistSeries,myROIi] = slimportROItseries(o,si,ROIi)

%This function load ROI Tseries data that have been loaded with
%loadROITSeries and then saved in a tSeries...file

%note: contains information for many sessions
%si : session index

%import roi data from input file path
[sPath,fn] = fileparts(o.sessPath{si});
if isfield(o,'SavedroiTSeriesfilename ')
    %import roi file from input path
    savedDir = o.SavedroiTSeriesfilename;    
    myROIi = importdata(savedDir);
    thistSeries = dir(savedDir);
    cd(o.sessPath{si})  
    return
    %mport roi data from default file path
else
    %move to default roi path
    savedDir = [o.stckPath 'slStckAnalyses/' o.myGroup '/' fn '/roiTseries/' o.myROIname{ROIi}];
end

%warning if roi tseries do not exist
if ~exist(savedDir)
    fprintf('%s \n','(slimportROItseries) Tseries do not exist for this ROI. Create them first by not removing "loadSavedROI" from varargin')
    fprintf('%s \n','(slimportROItseries) or the directory might be wrong ...')
    keyboard    
end
    cd(savedDir)
%check that there only is one file
f = dir([savedDir '/tSeries*']);

%get most recent file
lastfile = slGetMostRecentFileIndir;

%move othe file to archive
if length([f.isdir])>1;
    %move older files to archive
    mkdir archive    
    for i = 1 : length(f)
        if ~strcmp(f(i).name,lastfile)
            movefile(f(i).name,[savedDir '/archive'])
        end
    end
    fprintf('%s \n','(slImportROItseries) Getting most recent roi tseries file.')
end

%file
thistSeries = dir([savedDir '/' lastfile]);

%import roi data
myROIi = importdata([savedDir '/' thistSeries.name]);
cd(o.sessPath{si})