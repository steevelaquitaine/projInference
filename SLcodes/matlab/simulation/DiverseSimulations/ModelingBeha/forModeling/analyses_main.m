%% Prior direction distribution biases perception as predicted by Bayesian inference
% INTRODUCTORY TEXT: 
    % Steeve laquitaine
    % started 12/09/07 - now
    % usage: this script summarizes all the analyses performed

%% figure 1.Expected results
%
%<<120816_steeve_JSPS_expResults_lowRes_lowsz.PNG>>
%
%figure 3. Expected results.
 
%% figure 2.Prior distributions
%
%
%*wide prior (std = 1000 degrees)*
%open('steeve_exp07_metho_Pstd0020_mean225_card270_coh008_024_1_dir15_t62perCoh_121130_his.fig')
%open('steeve_exp07_metho_Pstd0020_mean225_card270_coh008_024_1_dir15_t62perCoh_121130_pol.fig')
%%
%*Narrow prior (std = 10 degrees)*
%open('steeve_exp07_metho_Pstdinf_mean225_card270_coh008_024_1_dir15_t60perCoh_121130_his.fig')
%open('steeve_exp07_metho_Pstdinf_mean225_card270_coh008_024_1_dir15_t60perCoh_121130_pol.fig')
%

%% figure 3. *A strong prior biases estimation when motion is weak*
tic
datlist.subjects={'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08'};
fig.nm='steeve_descriptive_analyses';
[figOut,GenStat]=analyses(datlist,{'coh','Pstd','sample_dir'},fig);
toc



%% figure 4. *Modeling*
%2900 seconds
tic
datlist.subjects={'sub01'};
fig.nm='steeve_modeling';
[Sp,Pred,fPa,R2,udata,sdata,dataOutput,PredOutput,FAbp,FA,StdPa]=modelBayesInf4(datlist,...
    {'coh','Pstd','sample_dir'},...
    fig,...
    [],[]);
toc
%problem with subject 4