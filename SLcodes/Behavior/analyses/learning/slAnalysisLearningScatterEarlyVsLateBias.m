

%slAnalysisLearningScatterEarlyVsLateBias.m
%
%
%
% author: steeve laquitaine
%   date: 160511
% purpose: Scatter plot the slopes of estimates vs true motion direction for
%         early versus late trials pooled across runs.
%         N Early and N late trials are sorted for each prior and pooled
%         together 
%         Trial data are linearized as data + distances to true motion 
%         direction (now range between -180 deg and 540) then         
%         averaged for each direction, coherence and prior condition.
%         The means are then fitted with a line producing a slope estimate
%         for each condition.
%         
%
% usage:
%
%   d = SLMakedatabank({'sub01'},'dataPath','~/data/dataPsychophy/proj01_priorStrength/','experiment','vonMisesPrior');
%   [linfit_early,linfit_late,r,p] = slAnalysisLearningScatterEarlyVsLateBias(d,100)
%
%Inputs:
%
%     d: database from "SLmakedatabank.m"
%   per: number of trials in early and late period of each run
% 
%outputs:
%
%       linfit_early: 
%            linfit_early.slope_fitMean : linear fit slopes for each prior/ 
%                                         coherence condition
%           linfit_early.interc_fitMean : linear fit intercept
%                      linfit_early.coh : coherence conditions
%                    linfit_early.prior : prior conditions
%       linfit_late: same but for late trials
%                 r: correlation between early and late slopes
%                 p: p-value for null hypothesis that slopes are not
%                    correlated

function [linfit_early,linfit_late,r,p] = slAnalysisLearningScatterEarlyVsLateBias(d,per)

%id 100 first and 100 last trials for each run
run = cell2mat(d.run);
run_u = unique(run);
firstlastID = nan(size(run));
%==================================================
%id early and late trials with 1 and 2 for each run
%==================================================
for i = 1 : length(run_u)
    runi_trials = find(run==run_u(i));
    firstlastID(runi_trials(1:per)) = 1;
    firstlastID(runi_trials(end-per+1:end)) = 2;
    runi_trials = [];
end

%==================================================
%Remove all successive pairs of blocks with same
%priors as next block. Early bias might be caused
%by previous block instead of quick learning.
%==================================================
%first run and second of runs with same contiguous priors
rcont = slgetSameSuccessiveRuns(d.Pstd,run);

%Analyse
priors_u = sort(unique(d.Pstd),'descend');
coh_u = sort(unique(d.stimStrength),'descend');
npri = length(priors_u);
ncoh = length(coh_u);
%find unique prior and coherence conditions
cnt = 0;

%get predictor (indep.:estimates) and predicted (depdt:directions)
%[feat,est] = slLinearizeCircDataForScatterAndLinFit(d.stimFeatureDeg,d.estimatesDeg);
feat = d.stimFeatureDeg;
est = d.estimatesDeg;

%prior mean
pmean = unique(d.priormean);

%loop over priors and coherences
%and plot early and late estimates
%against directions
pos_early = 1:2:ncoh+2;
pos_late = pos_early+1;

for pi = 1 : npri
    for ci = 1 : ncoh
        cnt = cnt+1; cnt_st = num2str(cnt);
        %========================
        %select early trials data
        %========================
        %each task condition and successive runs with different priors
        feat_early = feat(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==1 & rcont==1);
        est_early = est(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==1 & rcont==1);
        figure(1)
        np_st = num2str(npri); nc_st = num2str(ncoh);
        hold on; hline(unique(pmean),':b')
        subp = ['subplot(' np_st ',' nc_st ',' cnt_st ')'];
        %linear fits
        [~,~,linfit_early(pi,ci)] = SLdrawCircStat(est_early,feat_early,[],[],...
            'mysubplots',subp,'mycolor',[.3 .3 .3]);
        
        %=======================
        %select late trials data
        %=======================
        %same as above for late
        feat_late = feat(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==2 & rcont==1);
        est_late = est(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==2 & rcont==1);
        %linear fits
        [~,~,linfit_late(pi,ci)] = SLdrawCircStat(est_late,feat_late,[],[],...
            'mysubplots',subp,'mycolor',[1 0 0]);                
    end
end

%record conditions
for pi = 1 : npri
    for ci = 1 : ncoh
        linfit_early(pi,ci).coh = coh_u(ci);
        linfit_early(pi,ci).prior = priors_u(pi);
    end
end

%scatter slopes of fit to raw data
%priors/coherecnes are colored/scaled
figure('color','w')
%sort everything by coherence for plot
cohs = [linfit_early.coh];
[~,ix] = sort(cohs,'ascend');
cohs = cohs(ix);
priors = [linfit_early.prior];
priors = priors(ix);
%set scatter plot prior colors and coherence scales
c = [[0.5 0 0;1 0 0;1 0.5 0;0.65 0.65 0]];
colors = nan(length(priors),3);
for i = 1 : npri
    p_i = find(priors==priors_u(i));
    colors(p_i,:) = repmat(c(i,:),length(p_i),1);
end
sc = [400 200 100];
for i = 1 : ncoh
    scales(cohs==coh_u(i)) = sc(i);
end
%scatter plot
slopes_early = [linfit_early.slope_fitMean];
slopes_late = [linfit_late.slope_fitMean];
scatter(slopes_early(ix),slopes_late(ix),scales,colors,'filled')
hold on; scatter(slopes_early(ix),slopes_late(ix),scales,'w')
hold on; plot([0:1],[0:1],'k--')
ylim([-0.1 2])
xlim([-0.1 2])
xlabel({'Slopes early (linear fit slopes)','dot scales reflect coherence; colors, priors'})
ylabel('Slopes late (linear fit slopes)')

%======
%Stats:
%======
%Are early and late slopes strongly correlated ?
[r,p] = corrcoef(slopes_early(ix),slopes_late(ix));
title(['Correlation: R=' num2str(round(r(2,1)*100)/100) '; p=' num2str(p(2,1))])








