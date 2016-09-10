
%author: steeve laquitaine


function [Pdata,xpdf,o] = slmakeDataDistOverSub(subjects,db,dataBins,varargin)

%get arguments
if any(strcmp(varargin,'vonMisesPrior'))
    priorShape = 'vonMisesPrior';
end

numSub = length(subjects);
subjectsNum = unique(db.subjects);

%data and predictions for each subject
for sub = 1 : numSub
    
    %this trial's subject
    TrialsThisSub = subjectsNum(sub)==db.subjects;
    
    %data and conditions
    data = round(db.estimatesDeg(TrialsThisSub));
    d    = db.stimFeatureDeg(TrialsThisSub);
    coh  = db.stimStrength(TrialsThisSub);
    pstd = db.Pstd(TrialsThisSub);
    priorModes = cell2mat(db.priormodes(TrialsThisSub,:));
    
    %data mean and std
    [meanDatatmp,stdDatatmp,dataCondtmp] = SLmakeDataMeanAndStd(data,...
        d,coh,pstd,priorModes,priorShape);
    dataCondtmp = sortrows(dataCondtmp,[1 2 3]);
    
    %subject's estimate distribution
    %             [pData,xpdf] = makeDataDist(data,d,coh,pstd,priorModes,cond,...
    %                 priorShape);
    [Pdata,xpdf,o] = slmakeDataDist(dataBins,db.estimatesDeg,db.stimFeatureDeg,db.stimStrength,db.Pstd,[db.priormodes{:}],'vonMisesPrior');
    
    %make sure predicted and data distributions are calculated on
    %the same space. Works for predicted estimate space = 1:1:360 deg.
    %adjust predicted estimate probability distribution from 1:1:360 to
    %0:10:360 by summing probabilities within consecutive bins of 10
    %deg (law of probabilities).
    commonSpace = xpdf;
    [~,bins]= histc(0:1:360,commonSpace);
    bins(end) = [];
    bins = bins';
    commonSpace = commonSpace(2:end);
    
    %match conditions, data and predictions across predictions and
    %data and across subjects
    numDataEst = size(Pdata,1);
    for j = 1 : size(o.uniqCond,1)
        
        %case cond (predictions) and dataCondtmp (data) are same as
        %expected
        %get current condition sync with a reference ordering of the
        %conditions
        if isequal(dataCondtmp,o.uniqCond)
            
            %get this condition
            [~,posThisCond] = ismember(o.uniqCond(j,:),...
                dataCondtmp,'rows');
            
            %case condition exists
            if ~isequal(posThisCond,0)
                
                %match data and predictions mean, std and distributions
                %between conditions and subjects to the order in
                %o.uniqCond
                o.meanData(j,sub) = meanDatatmp(posThisCond);
                o.stdData(j,sub)  = stdDatatmp(posThisCond);
                o.PdisData(:,j,sub) = Pdata(:,posThisCond);
            else
                %else NaN
                o.meanData(j,sub) = NaN;
                o.stdData(j,sub)  = NaN;
                o.PdisData(:,j,sub) = NaN(numDataEst,1);
            end
        else
            fprintf('%s \n','(slmakeDataDistOverSub) Conditions do not match.')
        end
    end
end

%stop or continue ?
YesNo = input('Do you want continue: y/n ?','s');
if strcmp(YesNo,'n'); keyboard; end

%---------------------
%average over subjects
%---------------------

%calculate data mean, std and distributions
for i = 1 : size(o.meanData,1)
    CircMeanOvSub = SLcircMeanStd(o.meanData(i,:)','polar');
    o.meanDataOvSub(i) = CircMeanOvSub.deg.mean;
end
o.meanDataOvSub = SLmakeColumn(o.meanDataOvSub);
o.stdDataOvSub  = nanmean(o.stdData,2);
o.meanDisDataOvSub = nanmean(o.PdisData,3);
end