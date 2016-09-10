
%slSaveClassiRes.m
%
%
% author: steeve laquitaine
%purpose: save the results of classifications in a structured path
%
%  usage : 


%save classification data
function [o,c,s,roii] = slSaveClassiRes(o,c,s,roii)

%convert cell class to string
o.savedClass = o.myClass;
while iscell(o.savedClass)
    o.savedClass = [o.savedClass{:};];
end
%if not testing save data in structured dir
if o.test~=1
    %backup
    if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
    cd(o.myROIname{roii})
    o.myAnalPath2{roii} = pwd;
    %archive existing file
    list = dir('*.mat');
    if ~isempty(list)
        if isempty(dir('Archive')); mkdir Archive; end
        for i=1:length(list); movefile(list(i).name,'Archive/'); end
    end
    save(['ClassifStckSess' o.savedClass '_' o.myROIname{roii} '.mat'],'o','c','s')
    %save chance accuracy results
    %in a different file
    chance.raw = c.myClasf.raw.fullSgShuf;
    chance.Zsc = c.myClasf.Zsc.fullSgShuf;
    save(['ClassifStckSessChance' o.savedClass '_' o.myROIname{roii} '.mat'],'chance')
    clear chance;
    cd ..
        
elseif o.test==1
    %otherwise save in test folder
    %backup
    if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
    cd(o.myROIname{roii})
    mkdir test; cd test    
    o.myAnalPath2{roii} = pwd;
    %archive existing file
    list = dir('*.mat');
    if ~isempty(list)
        if isempty(dir('Archive')); mkdir Archive; end
        for i=1:length(list); movefile(list(i).name,'Archive/'); end
    end
    save(['ClassifStckSess' o.savedClass '_' o.myROIname{roii} '.mat'],'o','c','s')
    %save chance accuracy results
    %in a different file
    chance.raw = c.myClasf.raw.fullSgShuf;
    chance.Zsc = c.myClasf.Zsc.fullSgShuf;   
    save(['ClassifStckSessChance' o.savedClass '_' o.myROIname{roii} '.mat'],'chance')
    clear chance;
end

%cleanup
if isfield(c.myClasf.raw,'fullSgShuf')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShuf');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShuf');
end
if isfield(c.myClasf.raw,'fullSgShufBal')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShufBal');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShufBal');
end
