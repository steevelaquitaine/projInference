
%INITPRIORS Summary of this function goes here
%   Detailed explanation goes here
% [121126] - Steeve

% INPUT
    % mean
    % variance
    % trialnum

% OUTPUT
    % vector of n trials directions: task.parameter.dir.series


% a good example for our purpose: to choose a good std, we have to trade
% between:
% A)"setting a very narrow distribution and get prior's maximum effect, but
% keeping the full range of sample directions (lower probability directions may
% not occur at all").

% B)"setting the distribution width so that we always get at least 1 occurrence
% of the lower probability directions hence maintaining the full range of sample
% directions. But 1 occurence is not sufficient for statistics over the
% lower probability directions" --> no conclusions are possible for lower
% proba.directions

% C) "setting the distribution width so that we always get at least 5
% occurences of the lower proba.directions". The effect of the prior is reduced but
% we maintain a same full range of sample directions over conditions
% (cannot counfound). Moreover non parametric statistics becomes possible for lower probability
% directions.

% --> I chose A in this set of experiments to first check if we have any
% effect of the prior at all.

% plot the arrows on a unit circle
% modified from http://stackoverflow.com/questions/1803043/how-do-i-display-an-arrow-positioned-at-a-specific-angle-in-matlab

%--------------------------------------------------------------------------
% Updates
%--------------------------------------------------------------------------
% [121201]
% Current prior generation method is suceptible to noise due to random sampling
% e.g., distribution not peak at the mean or shape asymetry around the
% mean. I changed method. Now the number of repetition of each sample 
% is approximated from a continuous distribution. This method ensures a
% perfectly symetric gaussian distribution peaking at the mean that is not
% susceptible to noise in random sampling. 

%--------------------------------------------------------------------------
% Inputs
%--------------------------------------------------------------------------
% [121130]
% usage:
% You have to input the following parameters
% % Weak prior
% [task] = initPriors('steeve_exp07_metho_Pstd1000_mean225_card270_coh008_024_1_dir15_t60perCoh_121130',225,1000,60,155:10:295)
% note for a uniform distribution use:
% [task] = initPriors('steeve_exp07_metho_Pstdinf_mean225_card270_coh008_024_1_dir15_t60perCoh_121130',225,inf,60,155:10:295)
% % or Strong prior
% [task] = initPriors('steeve_exp07_metho_Pstd0020_mean225_card270_coh008_024_1_dir15_t60perCoh_121130',225,20,60,155:10:295)



% [121226] 
% % Weak prior (uniform)
% [task] = initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,inf,[.12 .35 1],[105 75 30],155:10:295)

% % Strong prior (gaussian)
% [task] = initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,20,[.12 .35 1],[105 75 30],155:10:295)



% [130202] 
% % No prior (uniform, all circle)
% [task] = initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,inf,[.12 .35 1],[105 75 30],0:10:350)

% % weak prior (slight bias)
% [task] = initRun('steeve_exp08_metho_Pstd100_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,100,[.12 .35 1],[105 75 30],155:10:295)

% % Intermediate weak prior (higher bias)
% [task] = initRun('steeve_exp08_metho_Pstd050_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,50,[.12 .35 1],[105 75 30],155:10:295)

% % clear prior (clear bias)
% [task] = initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,20,[.12 .35 1],[105 75 30],155:10:295)

% % Strong prior (strong bias)
% [task] = initRun('steeve_exp08_metho_Pstd14_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,14,[.12 .35 1],[105 75 30],155:10:295)


%% current test

% [130204]

% % Strong prior
% [task] = initRun('steeve_exp08_metho_Pstd14_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, 15,[.08 .16 32],[105 75 30],5:10:355)

% % Less strong prior
% [task] = initRun('steeve_exp08_metho_Pstd20_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, 20,[.08 .16 32],[105 75 30],5:10:355)

% % Weaker prior
% [task] = initRun('steeve_exp08_metho_Pstd40_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, 40,[.08 .16 32],[105 75 30],5:10:355)

% % No prior
% [task] = initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh008_016_032_dir36_t144_108_36perCoh_130205',225,367,[.08 .16 .32],[144 108 36],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh008_016_032_dir36_t108_072_36perCoh_130205',225,367,[.08 .16 .32],[108  72 36],5:10:355)



% [130205]

% % Strong prior
% [task] = initRun('steeve_exp08_metho_Pstde2_50_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, exp(1.25*2),[.08 .16 .32],[108 72 36],5:10:355)

% % Less strong prior
% [task] = initRun('steeve_exp08_metho_Pstde3_75_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, exp(1.25*3),[.08 .16 .32],[108 72 36],5:10:355)

% % Weaker prior
% [task] = initRun('steeve_exp08_metho_Pstde5_00_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, exp(1.25*4),[.08 .16 .32],[108 72 36],5:10:355)

% % No prior
% [task] = initRun('steeve_exp08_metho_Pstde6_25_mean225_card180270_coh008_016_032_dir36_t144_108_36perCoh_130205',225, exp(1.25*5),[.08 .16 .32],[144 108 36],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstde6_25_mean225_card180270_coh008_016_032_dir36_t108_072_36perCoh_130205',225, exp(1.25*5),[.08 .16 .32],[108  72 36],5:10:355)




% [130210]
% % Strong prior (even lower coherence)
% [task] = initRun('steeve_exp08_metho_Pstd012_mean225_card180270_coh004_012_028_dir36_t105_75_30perCoh_130209'  ,225, exp(1.25*2),[.04 .12 .28],[108 72 36],5:10:355)

% % Less strong prior
% [task] = initRun('steeve_exp08_metho_Pstd043_mean225_card180270_coh004_012_028_dir36_t105_75_30perCoh_130209'  ,225, exp(1.25*3),[.04 .12 .28],[108 72 36],5:10:355)

% % Weaker prior
% [task] = initRun('steeve_exp08_metho_Pstd149_mean225_card180270_coh004_012_028_dir36_t105_75_30perCoh_130209'  ,225, exp(1.25*4),[.04 .12 .28],[108 72 36],5:10:355)

% % No prior
% [task] = initRun('steeve_exp08_metho_Pstd518_mean225_card180270_coh004_012_028_dir36_t144_108_36perCoh_130209',225, exp(1.25*5),[.04 .12 .28],[144 108 36],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstd518_mean225_card180270_coh004_012_028_dir36_t108_072_36perCoh_130209',225, exp(1.25*5),[.04 .12 .28],[108  72 36],5:10:355)



% [130215]
% New priors (same coherence as previous)
% [task] = initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh004_012_028_dir36_t107_075_033perCoh_130215'  ,225, exp(1.25*2.4),[.04 .12 .28],[107 75 33],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh008_016_032_dir36_t107_075_033perCoh_130215'  ,225, exp(1.25*2.4),[.08 .16 .32],[107 75 33],5:10:355)


% [task] = initRun('steeve_exp08_metho_Pstd016_mean225_card180270_coh004_012_028_dir36_t107_075_032perCoh_130215'  ,225, 16,[.04 .12 .28],[107 75 32],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstd016_mean225_card180270_coh008_016_032_dir36_t107_075_032perCoh_130215'  ,225, 16,[.08 .16 .32],[107 75 32],5:10:355)


% [task] = initRun('steeve_exp08_metho_Pstd018_mean225_card180270_coh004_012_028_dir36_t107_075_032perCoh_130215'  ,225, 18,[.04 .12 .28],[108 73 33],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstd018_mean225_card180270_coh008_016_032_dir36_t107_075_032perCoh_130215'  ,225, 18,[.08 .16 .32],[108 73 33],5:10:355)


% [task] = initRun('steeve_exp08_metho_Pstd030_mean225_card180270_coh004_012_028_dir36_t107_075_032perCoh_130215'  ,225, 30,[.04 .12 .28],[108 73 33],5:10:355)
% [task] = initRun('steeve_exp08_metho_Pstd030_mean225_card180270_coh008_016_032_dir36_t107_075_032perCoh_130215'  ,225, 30,[.08 .16 .32],[108 73 33],5:10:355)


% [130216]
% [task] = initRun('steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_073_033perCoh_130216', 225, 10,[.06 .12 .24],[107 73 33],5:10:355)
% [task] = initRun('steeve_exp12_metho_Pstd020_mean225_coh006012024_dir36_t107_075_033perCoh_130216', 225, 20,[.06 .12 .24],[107 75 33],5:10:355)
% [task] = initRun('steeve_exp12_metho_Pstd040_mean225_coh006012024_dir36_t107_075_033perCoh_130216', 225, 40,[.06 .12 .24],[101 75 31],5:10:355)
% [task] = initRun('steeve_exp12_metho_Pstd080_mean225_coh006012024_dir36_t107_075_033perCoh_130216', 225, 80,[.06 .12 .24],[106 75 34],5:10:355)



% [130217]
% [task] = initRun('steeve_exp12_metho_Pstd080_mean225_coh008_dir36_t107_075_033perCoh_130217', 225, 80,[.08],[106],5:10:355)

%-------------------------------------------------------------------------
%%% Design a run
%-------------------------------------------------------------------------
function [task] = initRun(Prname,Prmean,Prstd,coh,Prnumtrials_perCoh,PrSp)
% Initialize figure
fig.hdle = figure(1);
fig.name = Prname; 

%%% Prnumtrials_perCoh: is a vector indicating the number of trials of
% different conditions (e.g., coherences).
Prior.parameter.dir.trialnum = [];
Prior.parameter.dir.series = [];
Priortmp = {nan};
Prior.parameter.dir.count = [];
Prior.parameter.dir.coh = [];
for i = 1:numel(Prnumtrials_perCoh)
    % Combine conditions of Priors (e.g., factor 1) and coherences (e.g., factor 2).
    [Priortmp{i}] = initPriors(Prname,Prmean,Prstd,Prnumtrials_perCoh(i),PrSp);
    % Initialize task
    task = Priortmp{1};
    % Pool directions (e.g., factor 3) across coherences.
    Prior.parameter.dir.trialnum = [Prior.parameter.dir.trialnum Priortmp{i}.parameter.dir.trialnum];
    Prior.parameter.dir.series =   [Prior.parameter.dir.series   Priortmp{i}.parameter.dir.series];
    Prior.parameter.dir.count =    [Prior.parameter.dir.count;   Priortmp{i}.parameter.dir.count];
    % Combine directions and coherences.
    %Prior.parameter.dir.coh=[Prior.parameter.dir.coh repmat(coh(i),1,Prnumtrials_perCoh(i))];
    Prior.parameter.dir.coh=[Prior.parameter.dir.coh repmat(coh(i),1, sum(Prior.parameter.dir.count(i,:),2))];
end
% Update task
task.parameter.dir.trialnum = sum(Prior.parameter.dir.trialnum,2);
task.parameter.dir.series   = Prior.parameter.dir.series;
task.parameter.dir.count    = sum(Prior.parameter.dir.count,1);
task.parameter.dir.coh      = Prior.parameter.dir.coh;

% Draw distributions as polar plots
[axr] = drawrosedir(task,fig);
% Draw distributions as histograms
[axh] = drawhistdir(task,fig);
% Backup figure
autobackup(fig.hdle,fig.name,'.fig')
% Backup parameters
autobackup(task,fig.name,'.mat')

%-------------------------------------------------------------------------
%%% Design Prior
%-------------------------------------------------------------------------
function [task] = initPriors(Prname,Prmean,Prstd,Prnumtrials,PrSp)
% fig.name=Prname; 
%mean of the distribution (stay the same).
task.parameter.dir.mean   = Prmean;   
%std of the distribution (stay the same).  
task.parameter.dir.std    = Prstd;  %std % low:10 (instead of 80); high:1000
%nb of trials for each coherence in a sub-block
task.parameter.dir.trialnum = Prnumtrials;%65; 
% you can choose to input or not a set of samples (if not, a default sample is chosen based on the sample size) 
% task.parameter.dir.sample.degree=[90 100 110 120 130 140 150 160 170 180]; % set of sample directions (the artifacted directions in experiment 1 and 2) 
task.parameter.dir.sample.degree = PrSp; % extend to include 90 degrees.
% task.parameter.dir.sample.degree=[40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220]; % extend to include cardinal directions
%nb of sample directions
task.parameter.dir.sampsiz  = numel(task.parameter.dir.sample.degree);
%nb of occurence of each displayed directions
task.parameter.dir.count = [];

% if no series of directions has been input
if ~isfield(task.parameter.dir,{'series'})
    
    % if a set of sample direction has not been input either
    if ~isfield(task.parameter.dir,{'sample'})    
        % create a default set of sample directions and a default series of
        % directions to be displayed on screen
        [task]=initSampleDir(task); % create a set of sample directions
        
        % Generate a discrete gaussian prior by random sampling (noisy).
        % [task]=initPriorDis_disc(task);
        % Generate a perfectly symetric gaussian prior
        [task]=initPriorDis_cont(task); % create a series of n trials directions
        disp(' "new directions have been initialized and drawn based on default parameters" ')
    else
        % if a set of displayed directions has been input.
        % Update the sample size
        task.parameter.dir.sampsiz = numel(task.parameter.dir.sample.degree);
        
        % Generate a discrete gaussian prior by random sampling (noisy).
        % [task]=initPriorDis_disc(task);
        % Generate a perfectly symetric gaussian prior
        [task]=initPriorDis_cont(task);
        disp(' "a series of n trials directions have been produced based on your input set of sample directions" ')
    end
    % Use input parameters.
else
    disp('the loaded directions statistics have been drawn')
end

%-------------------------------------------------------------------------
%%% Set sample states (e.g., 19 directions separated by 10 degrees)
%-------------------------------------------------------------------------
function [task] = initSampleDir(task)

% initialize sample directions
for i=1:task.parameter.dir.sampsiz
    
    % calculate the directions(%Angle of arrow in degree, from x-axis)
    task.parameter.dir.sample.degree(i)= 360/task.parameter.dir.sampsiz*i;         
    
    % Rotate directions away from the cardinals.
%     rotation = 15; % in degree
    % or not
    rotation = 0; % in degree
    task.parameter.dir.sample.degree(i)= task.parameter.dir.sample.degree(i) + rotation;

    % adjust the radian values according to the quadrant of the trigo.circle
    if task.parameter.dir.sample.degree(i)>360
       task.parameter.dir.sample.degree(i)= task.parameter.dir.sample.degree(i)-360;
    end
end
% sort the directions (adjusting the radian values creates an unsorted series of directions)
task.parameter.dir.sample.degree=sort(task.parameter.dir.sample.degree);

%-------------------------------------------------------------------------
%%% Produce a discrete gaussian prior 
%-------------------------------------------------------------------------
function [task] = initPriorDis_disc(task)
%Output
task.parameter.dir.series = [];
task.parameter.dir.count  = [];

% Generate a discrete gaussian prior (noisy due to sampling)
if task.parameter.dir.std~=inf
    % generate a time series of normally distributed directions.
    % We want to manipulate the distributions over the same range of
    % directions. Hence, while manipulating sigma we should always check that
    % all the sample directions of our range appears at least once!
    % task.parameter.dir.mean=227; %mean of the dist.
    % sigma=60; %std of the dist
    numerator=nan;
    denominator=task.parameter.dir.std*sqrt(2*pi);
    for i=1:task.parameter.dir.sampsiz;
        numerator(i)=exp((-(task.parameter.dir.sample.degree(i)-task.parameter.dir.mean)^2)/(2*task.parameter.dir.std^2));
    end
    f=numerator/denominator; %gaussian proba.dist.function.
    
    % fix a random number generator
    seed=1; randn('state',seed)
    
    % let's say we draw 'trialnum' random values from the set of sample directions
    task.parameter.dir.series = randsample(task.parameter.dir.sample.degree,task.parameter.dir.trialnum,true,f);
    
    % draw the discrete distribution of the directions we have generated
    value_dir=1:task.parameter.dir.sample.degree(end)';
    count_dir=histc(task.parameter.dir.series,[1:max(task.parameter.dir.sample.degree)])';
    for i =1:task.parameter.dir.sampsiz
        % find the directions we are interested in in the table
        theta_row=value_dir==task.parameter.dir.sample.degree(i);
        task.parameter.dir.count(i)=count_dir(theta_row);
    end
end
% Or generate a uniform prior
if task.parameter.dir.std==inf
    % check if "task.parameter.dir.trialnum" is a multiple of
    % "task.parameter.dir.sampsiz".
    if rem(task.parameter.dir.trialnum,task.parameter.dir.sampsiz)==0
        disp(['--- A UNIFORM distribution is being drawn ----'])
        task.parameter.dir.count=repmat(task.parameter.dir.trialnum/task.parameter.dir.sampsiz,...
            1,task.parameter.dir.sampsiz);
        task.parameter.dir.series=repmat(task.parameter.dir.sample.degree,task.parameter.dir.count(1),1);
        task.parameter.dir.series= task.parameter.dir.series(:);
        task.parameter.dir.series= task.parameter.dir.series';
    else
        disp(['"task.parameter.dir.trialnum" must be a multiple of "task.parameter.dir.sampsiz"'])
        return
    end
end
 
%-------------------------------------------------------------------------
%%% Produce a continuous gaussian prior
%-------------------------------------------------------------------------
function [task] = initPriorDis_cont(task)
%Output
task.parameter.dir.series = [];
task.parameter.dir.count  = [];

% Generate a continuous gaussian prior (reduced noise)
if task.parameter.dir.std~=inf
    % Prior's shape
    PriorShape = gauss_distribution(task.parameter.dir.sample.degree,task.parameter.dir.mean,task.parameter.dir.std);
    % Set the number of trials we need.
    f2 = PriorShape * task.parameter.dir.trialnum ;
    % plot(x,f2,'-o')
    for i = 1:numel(task.parameter.dir.sample.degree)
        % Repeat samples according to a gaussian distribution
        f2repet = [];
        f2repet = repmat(task.parameter.dir.sample.degree(i),round(f2(i)),1);
        % Store the repeated samples
        task.parameter.dir.series = [task.parameter.dir.series; f2repet];
    end
    task.parameter.dir.series=task.parameter.dir.series';
    %display data
    % hist(f2repetbkp,x); 
end
% Draw a discrete distribution of the directions 
for i =1:task.parameter.dir.sampsiz
    % Count sample directions
    task.parameter.dir.count(i)=numel(find(task.parameter.dir.series==task.parameter.dir.sample.degree(i)));
end
    
% Or generate a uniform prior
if task.parameter.dir.std==inf
    % check if "task.parameter.dir.trialnum" is a multiple of
    % "task.parameter.dir.sampsiz".
    if rem(task.parameter.dir.trialnum,task.parameter.dir.sampsiz)==0
        disp(['--- A UNIFORM distribution is being drawn ----'])
        task.parameter.dir.count=repmat(task.parameter.dir.trialnum/task.parameter.dir.sampsiz,...
            1,task.parameter.dir.sampsiz);
        task.parameter.dir.series=repmat(task.parameter.dir.sample.degree,task.parameter.dir.count(1),1);
        task.parameter.dir.series= task.parameter.dir.series(:);
        task.parameter.dir.series= task.parameter.dir.series';
    else
        disp(['"task.parameter.dir.trialnum" must be a multiple of "task.parameter.dir.sampsiz"'])
        return
    end
end

%-------------------------------------------------------------------------
%%% Draw distributions
%-------------------------------------------------------------------------
% Polar
function [axr] = drawrosedir(task,fig)
% Draw each sample direction
% axr.hdle = figure(2);
axr.hdle = subplot(1,2,1);
axr.name = strcat(fig.name,'_pol'); % tag the name with polar
title('Displayed directions (in degree)','fontsize',18)

% Initialize arrows for displayed directions
x = 0;                          %# X coordinate of arrow start
y = 0;                          %# Y coordinate of arrow start
% Set L at 90 to have identical scales for both histograms.
L = 65; %max(task.parameter.dir.count);                   %# Length of arrow
arrowsz = 5;
arrowidth = 0;
textdisp.y = 0;
textdisp.x = 0.1 * L;

% Draw an angle histogram (rose plot): Angle histogram is nice because the
% bias in the direction statistics is clearly represented. If we were to
% plot subject's produced directions, we would expect a deviation of all
% directions toward the more probable direction (attractor effect) when the distribution
% is biased and a reduction in this effect when the uncertainty in the prior increases.

% Convert the sample directions from degree to radians
task.parameter.dir.sample.rad = task.parameter.dir.sample.degree * pi / 180;

% Scale radius of polar plot
polar(0,L); 
hold on;

% Draw polar plots
% p=polar([task.parameter.dir.sample.rad task.parameter.dir.sample.rad(1)],...
%     [task.parameter.dir.count task.parameter.dir.count(1)]);
p = polar(task.parameter.dir.sample.rad, task.parameter.dir.count);
set(p,'LineWidth',2,'color','r')

% Remove unnecessary lines and text.
% Find the lines in the polar plot
h = findall(gcf,'type','line');
% Remove the handle for the polar plot line from the array
h(h == p) = [];
% Delete other lines
delete(h);
% Find all text objects in the polar plot
t = findall(gcf,'type','text');
% Delete text objects
delete(t);

% Draw the arrows
for i = 1:task.parameter.dir.sampsiz
    % calculate the coordinates of the arrows
    xEnd(i) = L*cos(task.parameter.dir.sample.degree(i)*pi/180);          %# X coordinate of arrow end
    yEnd(i) = L*sin(task.parameter.dir.sample.degree(i)*pi/180);          %# X coordinate of arrow end
    points = linspace(0,task.parameter.dir.sample.degree(i));     %# 100 points from 0 to theta
    % Draw arrows
    hold on;
    axis equal;
    arrow([0 0] , [xEnd(i) yEnd(i)] , arrowsz, [], [], arrowidth);       %# Plot arrow
    plot(0,0,'o','Markersize',10,'MarkerEdgeColor','k','MarkerFaceColor','w');  %# Plot point
    
end
% Annotate displayed directions
for i = 1:task.parameter.dir.sampsiz
    if xEnd(i)<0;
        dispgain.x=1.5;
    else
        dispgain.x=0;
    end
    text(xEnd(i)+sign(xEnd(i))*dispgain.x*textdisp.x,yEnd(i)-textdisp.y,num2str(task.parameter.dir.sample.degree(i)),'fontsize',16)
end
% Histogram
function [axh] = drawhistdir(task,fig)
% Initialize figure
% axh.hdle = figure(3);
axh.hdle = subplot(1,2,2);
axh.name = strcat(fig.name,'_his'); %tag the name with histogram
% Drawing settings
myBarWidth = 1;
% Draw data
bar(task.parameter.dir.sample.degree,task.parameter.dir.count,...
    'FaceColor',[1 0 0]   ,...
    'BarWidth' ,myBarWidth,...
    'EdgeColor',[1 1 1]);
hold on; 
% Replicate the first direction at the end of x-axis
bar_xend=2*task.parameter.dir.sample.degree(end)-task.parameter.dir.sample.degree(end-1);
bar(bar_xend,task.parameter.dir.count(1),...
    'FaceColor',[1 0 0],...
    'BarWidth' ,10     ,...
    'EdgeColor',[1 1 1]);
% Set x-axis unit
xunit=1:1:task.parameter.dir.sampsiz;
% Positon x-axis labels
set(gca,'fontsize',7,...
    'xtick',[task.parameter.dir.sample.degree(xunit) bar_xend],...
    'xticklabel',[task.parameter.dir.sample.degree(xunit) task.parameter.dir.sample.degree(1)])
xlabel('Displayed directions (in degree)')
% Legend y-axis
ylabel('Count')
% Title the histogram
title('Prior','fontsize',18)
xlim([min(task.parameter.dir.sample.degree)-10 bar_xend+10])
ylim([0 65]);

%-------------------------------------------------------------------------
%%% Generate a gaussian distribution
%-------------------------------------------------------------------------
function functionf = gauss_distribution(x, mu, s)
% x: predicted variable
% mu: mean
% s: standard deviation
% height: it is the height of the gaussian peak
% p1 = -.5 * ((x - mu)/s) .^ 2;
% set the shaping parameters of the gaussian
height = .5;

p1 = -height * ((x - mu)/s) .^ 2;
p2 = (s * sqrt(2*pi));
f = exp(p1) ./ p2;
% Convert from arbitrary unit to probabilities.
functionf=f/sum(f);

%-------------------------------------------------------------------------
%%% Backup figures and parameters
%-------------------------------------------------------------------------
function autobackup(task,filename,filetype)
% input:
% structure called 'fig' with 2 fields:
% - 'name': e.g., '120910_Steeve_exp02_resul_run8_var1000_mean195_coh024_t250_sess2'
% - 'hdle': the handle of the figure you want to backup

% check if the file name exists
r=dir;
clear i
nametocheck=strcat(filename,filetype);
for i=1:length(r); scanres(i) = strcmp(r(i).name,nametocheck); end

% if the file name exists, increment
i=0;
if ~isempty(find(scanres==1))
    while ~isempty(find(scanres==1))
        i=i+1;
        % increment the name
        filename=strcat(filename,'_0',num2str(i));
        nametocheck=strcat(filename,filetype);
        % check if the name exists already
        for j=1:length(r); scanres(j) = strcmp(r(j).name,nametocheck); end
    end
    errorms=[' "This filename exists already. The new filename is "',filename,'" '];
    disp(errorms)
end
if strcmp(filetype,'.fig')
    saveas(task,filename,'fig'); %name is a string
elseif strcmp(filetype,'.mat')
    save(filename,inputname(1)); %name is a string
end


