

% author: steeve
%   date: 150525
%purpose: show fit results in a table
%
%  usage: SLfitTable('AIC')
%
%
%varargin
%--------
%
%       'AIC': for maximum likelihood fit data.
%
%
%In progress:
% - automatically load data for each model, subject and store them in a
%   the table

function SLfitTable(varargin)

%status
fprintf('%s \n','(SLfitBayesianModel) Please set a data path ...')
dataPath = uigetdir(cd,'Pickup your project e.g., /dataPsychophy/Exp01...');
cd(dataPath)

%move to fit data directory
cd modelfit

%select data type ('AIC', R-squared fit, etc...)
if sum(strcmp(varargin,'AIC'))==1
    
    %check folder exist
    if SLexistFolder('AIC')==1 
        cd AIC
    else
        fprintf('(SLfitTable) I am sorry your folder "',folder,'"does not exist. \n')
    end
end


%get models
modelall = dir('model*');
for i=1:length(modelall)
    models{i}=modelall(i).name;
end
nModels = length(models);

%collect data per model
for i = 1 : nModels
    
    fprintf(['(SLfitTable) I am collecting "',models{i},'" fit data ... \n'])
    cd(models{i})
    
    %get subjects
    suball = dir('sub*');
    for j=1:length(suball)
        sub{i}=suball(j).name;
    end
    nSub = length(sub);
    
    %collect data per subj
    for j=1:nSub
        
        fprintf(['(SLfitTable) I am collecting "',sub{i},'" fit data ... \n'])
        cd(sub{i})
        
        dataThisSub = ls('datafit*')
        
    end

    
    
end

