

% slfmriStckGetInstancedb.m
%
%     author: steeve laquitaine
%    purpose: create instance database with data and variables stacked over sessions
%      usage: 
%
%initialize
%
%         o.roiname = roiname;
%         o.dataPath = ['~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/' o.roiname '/'];
%         sessPath{1} = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
%         sessPath{2} = '~/data/datafMRI/sltaskdotdirfmri05/s002520150923';
%         sessPath{3} = '~/data/datafMRI/sltaskdotdirfmri05/s002520150925';
%         o.filename = filename;

%stack database over sessions

function [o,d] = slfmriStckGetInstancedb(o)

for i = 1 : length(o.sessPath)
    clear d
    
    %get database this session
    o.thisSessPath = o.sessPath{i};    
    d = slfmriGetInstancedb({'myRandomCoh','myRandomDir','mySwitch'},o);
    save(['d' num2str(i)],'d')
    
    %stack over session
    dfieldnames = fieldnames(d);
    for f = 1 : length(dfieldnames)
        if i == 1
            dStck.(dfieldnames{f}) = [];
        end
        dStck.(dfieldnames{f}) = [dStck.(dfieldnames{f}); d.(dfieldnames{f})];
    end
end
d = dStck;
clear dStck; clear dfieldnames; clear i; clear f
save('d','d')