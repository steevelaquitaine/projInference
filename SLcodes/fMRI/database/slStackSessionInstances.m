
%slStackSessionInstances.m
%
%author: steeve laquitaine
%purpose: stack session instances or load stacked sessions instances from
%         analysis path
%
%  usage : 
%
%           
%           [o,si,s,c] = slStackSessionInstances(o,'CalcInstances')



%stack instances over sessions
function [o,si,s] = slStackSessionInstances(o,varargin)

%when we want to calculate instances from scratch
if slIsInput(varargin,'CalcInstances')  
    
    %get instances stack over session and save in directory
    [o,si,s] = slgetSessionStackedInstances(o);        
    
    %or load
elseif slIsInput(varargin,'loadSavedInstances')    
    %convert cell class to string
    o.savedClass = o.myClass;
    while iscell(o.savedClass)
        o.savedClass = [o.savedClass{:};];
    end
    cd([o.stckPath 'slStckAnalyses/' o.myGroup 'classif' o.savedClass [o.myAnalysis{:}]])    
    %check file exists
    l = dir(['ClassifStckSessData' o.myClass '_allRois_*']);
    if isempty(l)
        slPrintfStr('slStackSessionInstances',...
            [' ClassifStckSessData' o.savedClass '_allRois_..not found...'])
        dbstack
        keyboard
    else
        load(l.name)
    end
else
    slPrintfStr('slStackSessionInstances','Please input either ')
    slPrintfStr('slStackSessionInstances','- loadSavedInstances : to load existing instances')
    slPrintfStr('slStackSessionInstances','- CalcInstances : to calculate instances')
    dbstack
    keyboard
end