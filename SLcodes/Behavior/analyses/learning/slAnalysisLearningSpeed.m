
%slAnalysisLearningSpeed.m
%
%
% author: steeve laquitaine
%purpose: Analyse how prior-induced bias changes from early to late trials
%         within a block.
%         Remove all successive pairs of blocks with same priors as next
%         block early bias might be caused by previous block instead of
%         quick learning
%
%
%  usage:
%   e.g.,
%
%     d = SLMakedatabank({'sub01'},'dataPath','~/data/dataPsychophy/proj01_priorStrength/','experiment','vonMisesPrior');
%     [linfit_early,linfit_late,anov2_early,anov2_late] = slAnalysisLearningSpeed(d,100,1000,'bootstrapSlopeAndIntercept')
%
%
%Inputs:
%
%     d: database from "SLmakedatabank.m"
%   per: number of trials in early and late period of each run
%
%varargin:
%
%    'bootstrapSlopeAndIntercept': bootstrap slope and intercept 
%
%Outputs:
%
%   linfit_early: slope and intercept for the linear fit between estimates
%                 and directions of early trials
%   linfit_late : slope and intercept for the linear fit of late trials
%   anov2_early : anova 2 way to evaluate the effect of prior and coherence
%                 on estimate/displayed slope in early trials
%    anov2_late : same for late trials
%                 note: Anova2w is reported as F(df1,df2) = Fisher value
%                 df1 is degree of freedom between tested group
%                 df2 is total # of rows (observations) - # of groups (df error)
%
%
%Description:
%       Linear fit of raw data is not possible because data are circular
%       and certain data will look like outliers in the linear
%       space (e.g., 1 deg estimate for 360 deg motion direction)
%       but are not in the circular space. Data are thus transformed
%       to :
%           Transf. estimates = motion direction + shortest signed vector
%           distance between estimates and motion direction, then were
%           linearly fitted.



function [linfit_early,linfit_late,anov2_early,anov2_late] = slAnalysisLearningSpeed(d,per,nboot,varargin)

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
        %=======================
        %select early trials data
        %=======================
        %each task condition and successive runs with different priors
        feat_early = feat(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==1 & rcont==1);
        est_early = est(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==1 & rcont==1);
        figure(1)
        np_st = num2str(npri); nc_st = num2str(ncoh);
        hold on; hline(unique(pmean),':b')
        subp = ['subplot(' np_st ',' nc_st ',' cnt_st ')'];
        %plot
        [~,~,linfit_early(pi,ci)] = SLdrawCircStat(est_early,feat_early,[],[],...
            'mysubplots',subp,'mycolor',[.3 .3 .3],varargin{:});        
        
        %======
        %legend
        %======
        figure(1)
        if ci == 1
            ylabel({'Linearized estimates','(x+signed angle(y,x),','deg)'})
        end
        if pi == 1
            title([num2str(coh_u(ci)) '%coh'])
        end
        if pi == npri
            xlabel('Motion direction (deg)')
        end
        
        %=======================
        %select late trials data
        %=======================
        %same as above for late
        feat_late = feat(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==2 & rcont==1);
        est_late = est(d.Pstd==priors_u(pi) & d.stimStrength==coh_u(ci) & firstlastID==2 & rcont==1);
        %visualize with linear fit
        [~,~,linfit_late(pi,ci)] = SLdrawCircStat(est_late,feat_late,[],[],...
            'mysubplots',subp,'mycolor',[1 0 0],varargin{:});        
        hold on; vline(unique(pmean),':b')
        %======
        %legend
        %======
        if pi == 1
            title([num2str(coh_u(ci)) '%coh'])
        end
        if pi == npri
            xlabel('Motion direction (deg)')
        end
        if ci == 1
            ylabel('Estimate std (deg)')
        end
        
        %=======================================
        %Compare fit slope early vs. late trials
        %=======================================
        m_early = nanmean(linfit_early(pi,ci).slope_boot);
        m_late = nanmean(linfit_late(pi,ci).slope_boot);
        s_early = slMakeCI(linfit_early(pi,ci).slope_boot,.95);
        s_late = slMakeCI(linfit_late(pi,ci).slope_boot,.95);
        figure(3);set(gcf,'color','w');hold all;
        subplot(1,npri,pi)
        title(['Prior ' num2str(priors_u(pi)) ' deg'])
        hold all; bar(pos_early(ci),m_early,'facecolor',[.3 .3 .3],'edgecolor',[.3 .3 .3],'displayName','Early');
        bar(pos_late(ci),m_late,'r','edgecolor','r','displayName','Late');
        legend show
        myerrorbar(pos_early(ci),m_early,'yError',s_early,'MarkerSize=1','linewidth=1')
        myerrorbar(pos_late(ci),m_late,'yError',s_late,'MarkerSize=1','linewidth=1')
        ylim([0 1.4])
        xlim([0 pos_late(end)+0.5])
        hline(1,'--k')
        xlabel('Coherence (%)')
        if pi==1
            ylabel('Linear fit slope')
        end
        if ci==1
            set(gca,'xtick',pos_early+0.5,'xticklabel',coh_u)
        end
    end
end

%record conditions
for pi = 1 : npri
    for ci = 1 : ncoh
        linfit_early(pi,ci).coh = coh_u(ci);
        linfit_early(pi,ci).prior = priors_u(pi);
    end
end

%=====
%Stats
%=====
%Double effect of prior and coherence on slope(Anova 2-way) in early 
%and late trials
%Organize database for anova2 with nBoot (e.g.1000) reps
db_earlySlopes = [];
db_lateSlopes = [];
for pi = 1 : npri
    for ci = 1 : ncoh
        db_earlySlopes = [db_earlySlopes; linfit_early(pi,ci).slope_boot];       
        db_lateSlopes = [db_lateSlopes; linfit_late(pi,ci).slope_boot];       
    end
end
db_earlySlopes = reshape(db_earlySlopes,ncoh*nboot,npri);
db_lateSlopes = reshape(db_lateSlopes,ncoh*nboot,npri);
%run anova 2-way
[~,anov2_early] = anova2(db_earlySlopes,nboot);
[~,anov2_late] = anova2(db_lateSlopes,nboot);
