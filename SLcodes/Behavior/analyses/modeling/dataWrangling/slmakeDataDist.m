
%slmakeDataDist.m
%
%
%author: steeve laquitaine
%
%
%
%
% usage:
%
%               data : series of data (e.g., direction estimates)
%                  d : associated factor 1 (e.g., motion direction)
%                coh : associated factor 2 (e.g., coherence)
%               pstd : associated factor 3 (e.g., prior std)
%         priorModes : associated factor 4 (e.g., prior modes)
%               cond : matrix of Ncond by 3 factors unique conditions (1 to 3)
%         priorShape : 'vonMisesPrior' or 'bimodalPrior'
% 
%         [p,x] = slmakeDataDist(data,d,sStrg,pstd,priorModes,priorShape)
%
%outputs 
%
%                  p: NdataBins-1 by Ncond matrix of probabilities
%                  p: 1by Bins

function [p,xpdf,o] = slmakeDataDist(dataBins,data,d,StimStrength,pstd,priorModes,priorShape,varargin)

%get task conditions
[o,nCond] = slGetMotionTaskCond(d,StimStrength,pstd,priorModes,priorShape);
cond = o.uniqCond;
o.nCond = nCond;

%start from 0, end at 360 
%to get full range
motDir = dataBins;

%this stimStrength
p = nan(numel(motDir)-1,size(cond,1));

%case bimodal prior
if strcmp(priorShape,'bimodalPrior')
    priorCond = priorModes(:,2) - priorModes(:,1);
    %sanity check
    numPriorCond = size(SLuniqpair(priorModes),1);
    if numel(unique(priorCond)) ~= numPriorCond
        fprintf('%s \n', ['(makeDataDist) Something wrong with number of',...
            'prior conditions for bimodal prior'])
    end
end

%this condition
%not on cluster
if ~any(strcmp(varargin{:},'sherlock')) && ~slIsInput(varargin{:},'noInput')
    YesNo = input('(slmakeDataDist) Do you want to see sample size for each condition y/n ?','s');
elseif slIsInput(varargin{:},'sherlock') || slIsInput(varargin{:},'noInput')
    %on cluster
    YesNo = 'y';
end

for i = 1 : size(cond,1)
    
    %case von Mises prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')
        thisCon = pstd==cond(i,1) & StimStrength==cond(i,2) & d==cond(i,3);
        
        %case bimodal prior
        %--------------------
    elseif strcmp(priorShape,'bimodalPrior')
        thisCon = priorCond==cond(i,1) & StimStrength==cond(i,2) & d==cond(i,3);
        
    else
        fprintf('%s \n', '(makeDataDist) You need to input prior type')
    end
    
    %get data for this condition
    dtoHist = data(thisCon);
    
    %make probability distribution
    [~,xpdf,ctmp(:,i)] = makePdf(motDir,dtoHist,'raw');
    
    %adjust count
    c(:,i) = ctmp(1:end-1,i);
    
    %probability
    p(:,i) = c(:,i)/sum(c(:,i));
    %[p(:,i),xpdf,c(:,i)] = makePdf(motDir,dtoHist,'raw');
    
    %if show sample size
    if strcmp(YesNo,'y')
        fprintf('%s \n', ['(makeDataDist) ',num2str(numel(dtoHist)),' trials'])
    end
    
end
