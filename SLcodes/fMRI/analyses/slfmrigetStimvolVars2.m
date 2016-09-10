%slfmrigetStimvolVars2.m
%
%
% author: steeve laquitaine
%   date: 160112
%purpose: retrieve variables associated with instances using stimvols 
%         (from getStimvol)
%
%usage:
%
%e.g., 
%     
%     %set params
%     %o.sessPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
%     o.sessPath = '/Volumes/Transcend/data/sltaskdotdirfmri05/s02520150814/';
%     %o.myGroup  = 'Concatenation';
%     o.myGroup  = 'MotionComp';
%     o.myScan     = 1;
%     o.taskNum    = 2;
%     o.phaseNum   = 1;
%     o.segmentNum = 2;
%     variables = {'myRandomDir','myRandomCoh'};
%       
%     %get stimvols you want variables of
%     cd(o.sessPath)
%     v = mrLoadRet([]);
%     v = viewSet(v,'curGroup',o.myGroup,'curScan',o.myScan);;
%     stimvols = getStimvol(v,'myRandomCoh','taskNum',o.taskNum,'phaseNum',...
%           o.phaseNum,'segmentNum',o.segmentNum);
%       
%     %get variable values associated with stimvols
%     [vars,o,db] = slfmrigetStimvolVars2(stimvols,variables,o)
% 
%
%inputs: 
%     stimvols: input stimvols
%    variables: variable values we want for stimvols
%            o: fmri session parameters
%
%output:           
%    vars: variables associated with stimvols
%      db: sorted variables and stimvols
%       o: parameters

function [vars,o,db] = slfmrigetStimvolVars2(stimvols,variables,o)

%------- get session scan ------
%open session
cd(o.sessPath)
v = mrLoadRet([]);

%set group & scan
v = viewSet(v,'curGroup',o.myGroup);
v = viewSet(v,'curScan',o.myScan);

dbstack
keyboard

%get all stimvols and variables
[d,o.variables] = getStimvol(v,variables,'taskNum',o.taskNum,'phaseNum',...
    o.phaseNum,'segmentNum',o.segmentNum);







