

%slLoadSavedFitParams.m
%
%
% author: steeve laquitaine
%purpose: load saved fit parameters and data for the different models of 
%         motion direction estimation
%
%usage:
%
%        path = '/Volumes/DroboBKUP/data/dataPsychophy/proj01_priorStrength/modelfit/AIC/model_Bayes_MAP_FatTailPrior/';
%        [AICs,allVars,sub] = slLoadSavedFitParams('allVars',path)
%
%
%inputs:
%     'AIC'/'allVars' : output AIC computed from all the trials or
%                            outputs all the variables for later 
%                            calculations

function [paraVals,allVars,sub] = slLoadSavedFitParams(para,path)

%go to path
cd(path)
subFd = dir('sub*');
nsub = length(subFd);

%preallocate
paraVals = nan(1,nsub);
allVars = [];
sub = nan(1,nsub);

for i = 1 : nsub
    cd(subFd(i).name)
    fitfile = dir('data*');
    %load file
    if ~isempty(fitfile)
        ws = load(fitfile.name);
    else
        fprintf('%s \n',['file for sub' subFd(i).name 'not found'])
        keyboard        
    end
    %case 'AIC' get AIC
    if strcmp(para,'AIC')
        %find logl per trial
        if isfield(ws,'Logl_pertrialBestfit')            
            logl = sum(ws.Logl_pertrialBestfit);
            nfitp = sum(~isnan(ws.fitP.p));
        elseif isfield(ws,'logl_pertrialBestfit')
            logl = sum(ws.logl_pertrialBestfit);
            nfitp = sum(~isnan(ws.fitP.p));
        elseif isfield(ws,'o')
            logl = ws.o.logl;            
            nfitp = sum(~isnan(ws.o.fitP));
        end                    
        AIC = 2*(nfitp - logl);
        paraVals(i) = AIC;
        sub(i) = i;
    end
    
    %case all parameters for later calculations
    if strcmp(para,'allVars')
        %find logl per trial
        if isfield(ws,'Logl_pertrialBestfit')            
            allVars(i).loglpertrial = ws.Logl_pertrialBestfit;
            allVars(i).nfitp = sum(~isnan(ws.fitP.p));
            allVars(i).disp = ws.output.disp;
            allVars(i).coh = ws.output.coh;
            allVars(i).pstd = ws.output.pstd;
        elseif isfield(ws,'logl_pertrialBestfit')
            allVars(i).loglpertrial = ws.logl_pertrialBestfit;
            allVars(i).nfitp = sum(~isnan(ws.fitP.p));
            allVars(i).disp = ws.output.disp;
            allVars(i).coh = ws.output.coh;
            allVars(i).pstd = ws.output.pstd;
        elseif isfield(ws,'o')
            allVars(i).loglpertrial = ws.o.logl_pertrialBestfit;            
            allVars(i).nfitp = sum(~isnan(ws.o.fitP));
            allVars(i).disp = ws.o.disp;
            allVars(i).coh = ws.o.coh;
            allVars(i).pstd = ws.o.pstd;
        end                            
        sub(i) = i;        
    end
    
    %clearex('subFd','nsub','paraVals')
    clear ws; clear logl; clear nfitp; clear AIC
    cd .. 
end