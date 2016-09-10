
%slwrapperPlotPriorBiasLearning.m
%
%
%
% author: Steeve Laquitaine
%   date: 16/03/08
%purpose: analyse the effect of the previous prior block (first column) 
%         on the prior-induced bias learnt in the next prior block 
%         (subsequent column). Figure columns alternate as follows :
%         previous(80-40-20-10) - next (one prior) - previous....
%
%Description : motion direction estimation task with 4 priors over motion directions
%              prior strengths are 80,40,20,10 deg in 4 randomly
%              interleaved blocks
%              plot estimates against true directions for 100st and last trials of each
%              prior block.Plot this for each 12 pairs of block sequences. 
%              e.g., 80 -> 40; 80 --> 20; 80 --> 10 ; 80 -->80 etc...
%              note: coherence are mixed

%% init analysis parameters
t1 = tic;
prevP      = [80 40 20 10]; %previous block prior
nextP      = [80 40 20 10]; %next block
N          = 100;
subjects   = {'sub01','sub02','sub03','sub04','sub05',...
    'sub06','sub07','sub08','sub09','sub10','sub11','sub12'};
datapath   = '~/data/dataPsychophy/proj01_priorStrength/';
experiment = 'vonMisesPrior';

%% get data
d = SLMakedatabank(subjects,'dataPath',datapath,'experiment',experiment);

%% Get Nst and last trials each prior block
%All trials are used for the linear fit (across subjects).
%Estimates are never averaged over any condition.
figure('color','w')
dotcolorst = [.8 .8 .8];
dotcoloren = [0.6 .6 .6];
fitcolorst = [.5 .5 .5];
fitcoloren = [0 0 0];

dotcolorstN = [1 .8 .8];
dotcolorenN = [1 .6 .6];
fitcolorstN = [1 .5 .5];
fitcolorenN = [1 0 0];

%for fit
linecolor=[0.7 0.2 0,1 0.2 0,1 0.7 0,0.7 0.7 0];

ax = 0;
for prev = 1:length(prevP)
    for next = 1:length(nextP)
        clear o        
        o = slgetDataBlockSequence(d,prevP(prev),nextP(next),100);

        %plot
        %=============
        %current block
        %=============
        %previous block - first trials
        ax = ax + 1;
        subplot(length(prevP),2*length(nextP),ax); hold all
        plot(o.dispPrevBNst,o.estPrevBNst,'.','color',dotcolorst,'markersize',10)
        
        %prev block - last trials
        plot(o.dispPrevBNend,o.estPrevBNend,'o','color',dotcoloren,'markersize',10)
        plot(1:360,1:360,'--k')
        hline(225,'--b')
        xlim([1 360]); ylim([1 360])
        title(['Previous block' num2str(prevP(prev))])

        %fits
        sllinefit(o.dispPrevBNst,o.estPrevBNst,1:1:360,fitcolorst,'linewidth',6)
        sllinefit(o.dispPrevBNend,o.estPrevBNend,1:1:360,fitcoloren,'linewidth',6)

        if ax==1
            ylabel({['Estimates (grey(' num2str(N) 'first'] ; ['(black:' num2str(N) 'last']})
        end

        %==========
        %Next block
        %==========
        %next - st
        ax = ax + 1;
        subplot(length(prevP),2*length(nextP),ax); hold all
        plot(o.dispNextBNst,o.estNextBNst,'.','color',dotcolorstN,'markersize',10)
        
        %next - end
        plot(o.dispNextBNend,o.estNextBNend,'o','color',dotcolorenN,'markersize',10)
        plot(1:360,1:360,'--k')
        hline(225,'--b')
        xlim([1 360]); ylim([1 360])
        title(['Next block' num2str(nextP(next))])
        
        %fits
        sllinefit(o.dispNextBNst,o.estNextBNst,1:1:360,fitcolorstN,'linewidth',6)
        sllinefit(o.dispNextBNend,o.estNextBNend,1:1:360,fitcolorenN,'linewidth',6)
        
        if ax==length(nextP)*2*length(nextP)
            xlabel({['Motion directions (deg); note: estimates were re-expressed with thei angle values',...
                'closest to the motion direction'],['[grey/black]: early/late trials previous black,[pink/red]:early/late next block']})
        end   
    end
end

elapsed = toc(t1);
disp(['Elapsed: ' elapsed ' sec'])