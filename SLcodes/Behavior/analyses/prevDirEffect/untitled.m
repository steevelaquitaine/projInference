


[SigBias2priorAll, SigBias2prevdir, marginAll,ciAll,mAll,sAll] = slplotandStatsPrevDirvsPriorBias(subject,'~/data/dataPsychophy/proj01_priorStrength/')

function [SigBias2priorAll, SigBias2prevdir, marginAll,ciAll,mAll,sAll]  = slplotandStatsPrevDirvsPriorBias(subject,datapath)

%% polar plots with average for each subject with stats over all conditions
[fig,GenStat,meanDispDirInfo,meanEstimateInfo,meanPrevDirInfo,FA1,FA2] = analyses(subject,...
    {'StimStrength','Pstd','FeatureSample'},...
    'dataPath','~/data/dataPsychophy/proj01_priorStrength/',...
    'experiment','vonMisesPrior',...
    'slgetprevDirEffect');

% %% convert to distances between -180 and 180 degs to be able to do statistics
% %is average distance between estimate and displayed direction 
% %   est ccw to disp --> distance > 0 ; if 0 ccw upper CI (biased to prior)
% %   est cw to disp --> distance < 0 ; if 0 cw lower CI (biased to previous)
% for i = 1 : size(meanDispDirInfo,1)
%    for j = 1 : size(meanDispDirInfo,2)       
%        %distance estimate to current displayed direction
%        [margin(i,j),ci,m(i,j),s(i,j)] = slMakeCI(SLvectors2signedAngle(meanEstimateInfo{i,j}.deg.all,meanDispDirInfo{i,j}.deg.all,'polar'),.95); 
%        cilow(i,j) = ci(1);
%        ciup(i,j) = ci(2);       
%    end
% end
% 
% %find conditions where the bias toward the prior (distance >0) is 
% %significant (lower CI > 0)
% posSigBias2prior = find(m>0 & cilow > 0);
% ciSigBias2prior = [cilow(posSigBias2prior)  ciup(posSigBias2prior)];
% meanSigBias2prior = m(posSigBias2prior);
% 
% % find conditions where the bias toward the previous direction(distance<0) is 
% %significant (higher CI < 0)
% posSigBias2prevdir = find(m<0 & ciup < 0);
% ciSigBias2prevdir = [cilow(posSigBias2prevdir)  ciup(posSigBias2prevdir)];
% meanSigBias2prevdir = m(posSigBias2prevdir);

%across all 12 conditions
estimates = [];
displ = [];
for i = 1 : size(meanDispDirInfo,1)
    for j = 1 : size(meanDispDirInfo,2)
        estimates = [estimates; meanEstimateInfo{i,j}.deg.all];
        displ     = [displ  ; meanDispDirInfo{i,j}.deg.all];
    end
end
[marginAll,ciAll,mAll,sAll] = slMakeCI(SLvectors2signedAngle(estimates,displ,'polar'),.95);

%is the overall bias toward the prior (distance > 0) significant (lower CI
%> 0) ?
SigBias2priorAll = mAll>0 & ciAll(1) > 0;

%is the overall bias toward the prior (distance < 0) significant (upper CI
%< 0) ?
SigBias2prevdir = mAll<0 & ciAll(2)<0;

fprintf('\n %s %i \n','Overall bias significance toward the prior : ',SigBias2priorAll)
fprintf('%s %.02f %s %.02f %.02f %s \n','mean distance to current direction:',mAll,', 95% CI [',ciAll,']')
fprintf('%s \n','note that distance > 0 indicate bias toward the prior ')
