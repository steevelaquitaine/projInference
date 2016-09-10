
%SLinitRunExpBimodalPrior.m
%
% author: Steeve Laquitaine
%   date: 150109
%purpose: Create the parameters of a motion direction experiment where:
%           - Each trial motion direction is drawn from a bimodal
%           distribution (mixture of von Mises)
%           - Motion coherences are randomized across trials
%
%usage:
%
%       task = SLinitRunExpBimodalPrior('steeve_exp02_Bimodal_metho_Pmodes205_245_Pstd020_mean225_coh6_12_24_t107_75_33_dir36_130412', 225, [145 305], 8.4,[.06 .12 .24],[107 75 33],5:10:355)
%
%Example of three bimodal priors tested in our experiment
%
% k = 23.75
% task = SLinitRunExpBimodalPrior('steeve_exp_BimodalVM_metho_Pmodes145_305_Pstd12_mean225_coh6_12_24_t107_75_33_dir36_150117',225,[145 305],23.75,[.06 .12 .24],[107 75 33],5:10:355)
% task = SLinitRunExpBimodalPrior('steeve_exp_BimodalVM_metho_Pmodes165_285_Pstd12_mean225_coh6_12_24_t107_75_33_dir36_150117',225,[165 285],23.75,[.06 .12 .24],[107 75 33],5:10:355)
% task = SLinitRunExpBimodalPrior('steeve_exp_BimodalVM_metho_Pmodes185_265_Pstd12_mean225_coh6_12_24_t107_75_33_dir36_150117',225,[185 265],23.75,[.06 .12 .24],[107 75 33],5:10:355)
%
%
%SLcodes required:
%
%   SLmixtureVM
%   SLmakeDiscreteMixtureVM


function task = SLinitRunExpBimodalPrior(Prname, Prmean,Prmodes,Prstrength, coh, Prnumtrials_perCoh, PrSp)

%figure
fig.hdle = figure(1);
set(gcf,'color','w','position',[440 544 560 254])
fig.name = Prname; 

%number of trials for each experimental condition 
%(e.g., motion coherences).
Prior.parameter.dir.trialnum = [];
Prior.parameter.dir.series = [];
Priortmp = {nan};
Prior.parameter.dir.count = [];
Prior.parameter.dir.coh = [];

for i = 1 : numel(Prnumtrials_perCoh)
    
    %Combine conditions of Priors (e.g., factor 1) and coherences (e.g., factor 2).
    [Priortmp{i}] = initPriors(Prname, Prmean, Prmodes, Prstrength, Prnumtrials_perCoh(i), PrSp);
    
    %Initialize task
    task = Priortmp{1};
    
    %Pool directions (e.g., factor 3) across coherences.
    Prior.parameter.dir.trialnum = [Prior.parameter.dir.trialnum Priortmp{i}.parameter.dir.trialnum];
    Prior.parameter.dir.series =   [Prior.parameter.dir.series   Priortmp{i}.parameter.dir.series];
    Prior.parameter.dir.count =    [Prior.parameter.dir.count;   Priortmp{i}.parameter.dir.count];
    
    %Combine directions and coherences.
    %Prior.parameter.dir.coh=[Prior.parameter.dir.coh repmat(coh(i),1,Prnumtrials_perCoh(i))];
    Prior.parameter.dir.coh = [Prior.parameter.dir.coh repmat(coh(i),1, sum(Prior.parameter.dir.count(i,:),2))];
end

%task parameters
task.parameter.dir.trialnum = sum(Prior.parameter.dir.trialnum,2);
task.parameter.dir.series   = Prior.parameter.dir.series;
task.parameter.dir.count    = sum(Prior.parameter.dir.count,1);
task.parameter.dir.coh      = Prior.parameter.dir.coh;

%polar plots of distributions
drawrosedir(task,fig);

%distributions as histograms
drawhistdir(task,fig);


%Backup
autobackup(fig.hdle,fig.name,'.fig')
autobackup(task,fig.name,'.mat')


%Trial distribution of direction
function task = initPriors(Prname, Prmean, Prmodes, Prstrength, Prnumtrials, PrSp)

%parameters
task.parameter.dir.ParaFilename = Prname;
task.parameter.dir.mean = Prmean;   
task.parameter.dir.modes = Prmodes;   
task.parameter.dir.strength = Prstrength;%von mises strength 

%nb of trials for each coherence in a sub-block
task.parameter.dir.trialnum = Prnumtrials;%65; 

% you can choose to input or not a set of samples (if not, a default sample is chosen based on the sample size) 
task.parameter.dir.sample.degree = PrSp;

%nb of sample directions
task.parameter.dir.sampsiz = numel(PrSp);

%nb of occurence of each displayed directions
task.parameter.dir.count = [];

%if no series of directions has been input
if ~isfield(task.parameter.dir,{'series'})
    
    % if a set of sample direction has not been input either
    if ~isfield(task.parameter.dir,{'sample'})    
        
        %default sample directions default series of directions
        task = initSampleDir(task);
        
        %status
        sprintf('(SLinitRunExpBimodalPrior) "new directions have been created and drawn based on default parameters" ')
    
    else
        %if a set of displayed directions has been input update
        task.parameter.dir.sampsiz = numel(task.parameter.dir.sample.degree);
        
        %status
        sprintf('(SLinitRunExpBimodalPrior) "A series of n motion directions have been produced" ')
    end
    task = SLmakeDiscreteMixtureVM(Prstrength,...
            PrSp,...
            Prmodes,...
            Prnumtrials);
        
    %Use input parameters.
else
    sprintf('(SLinitRunExpBimodalPrior) The input directions statistics have been drawn')
end

%Direction space
function task = initSampleDir(task)

% initialize sample directions
for i = 1 : task.parameter.dir.sampsiz
    
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


%Case noisy distribution (from sampling)
function [task] = initPriorDis_disc(task)
%Output
task.parameter.dir.series = [];
task.parameter.dir.count  = [];

% Generate a discrete gaussian prior (noisy due to sampling)
if task.parameter.dir.strength~=0
    % generate a time series of normally distributed directions.
    % We want to manipulate the distributions over the same range of
    % directions. Hence, while manipulating sigma we should always check that
    % all the sample directions of our range appears at least once!
    % task.parameter.dir.mean=227; %mean of the dist.
    % sigma=60; %std of the dist
    numerator=nan;
    denominator=task.parameter.dir.strength*sqrt(2*pi);
    for i=1:task.parameter.dir.sampsiz;
        numerator(i)=exp((-(task.parameter.dir.sample.degree(i)-task.parameter.dir.mean)^2)/(2*task.parameter.dir.strength^2));
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
if task.parameter.dir.strength==0
    % check if "task.parameter.dir.trialnum" is a multiple of
    % "task.parameter.dir.sampsiz".
    if rem(task.parameter.dir.trialnum,task.parameter.dir.sampsiz)==0
        sprintf(['(SLinitRunExpBimodalPrior)--- A UNIFORM distribution is being drawn ----'])
        task.parameter.dir.count=repmat(task.parameter.dir.trialnum/task.parameter.dir.sampsiz,...
            1,task.parameter.dir.sampsiz);
        task.parameter.dir.series=repmat(task.parameter.dir.sample.degree,task.parameter.dir.count(1),1);
        task.parameter.dir.series= task.parameter.dir.series(:);
        task.parameter.dir.series= task.parameter.dir.series';
    else
        sprintf(['(SLinitRunExpBimodalPrior) "task.parameter.dir.trialnum" must be a multiple of "task.parameter.dir.sampsiz"'])
        return
    end
end



%Draw 
%----
%Polar
function [axr] = drawrosedir(task,fig)

%Draw each sample direction
axr.hdle = subplot(1,2,1);
axr.name = strcat(fig.name,'_pol'); % tag the name with polar
title('Displayed directions (in degree)','fontsize',18)

% Initialize arrows for displayed directions
x = 0; 
y = 0; 

%Set L at 90 to have identical scales for both histograms.
L = 65;
arrowsz = 2;
arrowidth = 0;

%Draw an angle histogram (rose plot): Angle histogram is nice because the
%bias in the direction statistics is clearly represented. If we were to
%plot subject's produced directions, we would expect a deviation of all
%directions toward the more probable direction (attractor effect) when the distribution
%is biased and a reduction in this effect when the uncertainty in the prior increases.

% Convert the sample directions from degree to radians
task.parameter.dir.sample.rad = task.parameter.dir.sample.degree * pi / 180;

%Scale radius
polar(0,L); 
hold on;

%polar plots
p = polar(task.parameter.dir.sample.rad, task.parameter.dir.count);

%set(p,'LineWidth',2,'color','r','markeredgecolor','w')
set(p,'LineWidth',2,'color','r')

%Remove unnecessary lines and text.
h = findall(gcf,'type','line');

%Remove the handle for the polar plot line from the array
h(h == p) = [];

%Delete other lines
delete(h);

%Find all text objects in the polar plot
t = findall(gcf,'type','text');

%Delete text objects
delete(t);

%Draw arrows
for i = 1 : task.parameter.dir.sampsiz
    
    %arrows coordinates
    xEnd(i) = L*cos(task.parameter.dir.sample.degree(i)*pi/180);          %# X coordinate of arrow end
    yEnd(i) = L*sin(task.parameter.dir.sample.degree(i)*pi/180);          %# X coordinate of arrow end
    points = linspace(0,task.parameter.dir.sample.degree(i));     %# 100 points from 0 to theta
    
    %Draw arrows
    hold on;
    axis equal;
    arrow([0 0] , [xEnd(i) yEnd(i)] , arrowsz, [], [], arrowidth);       %# Plot arrow
    plot(0,0,'o','Markersize',10,'MarkerEdgeColor','k',...
        'MarkerFaceColor','w'); 
    
end

%Annotate directions
for i = 1 : task.parameter.dir.sampsiz
    
    text(xEnd(i)*1.1, yEnd(i)*1.1,...
        num2str(task.parameter.dir.sample.degree(i)),...
        'HorizontalAlignment','center','fontsize',9)
end

xlim([-65 65])
ylim([-65 65])

%Histogram
function [axh] = drawhistdir(task,fig)

%figure
axh.hdle = subplot(1,2,2);
axh.name = strcat(fig.name,'_his');

%graphics
myBarWidth = 1;

%draw
bar(task.parameter.dir.sample.degree, task.parameter.dir.count,...
    'FaceColor',[1 0 0]   ,...
    'BarWidth' ,myBarWidth,...
    'EdgeColor',[1 1 1]);
hold on; 

%x-axis
xunit = 1 : 5 : task.parameter.dir.sampsiz;
set(gca,'fontsize',7,...
    'xtick',task.parameter.dir.sample.degree(xunit),...
    'xticklabel',task.parameter.dir.sample.degree(xunit),...
    'fontsize',10)
xlabel('Displayed directions (in degree)','fontsize',12)

%y-axis
ylabel('Count','fontsize',12)

%Title
title('Prior','fontsize',18)
xlim([0 360])
ylim([0 65]);
box off



%Backup figures and parameters
%------------------------------
function autobackup(task,filename,filetype)
% input:
% structure called 'fig' with 2 fields:
% - 'name': e.g., '120910_Steeve_exp02_resul_run8_var1000_mean195_coh024_t250_sess2'
% - 'hdle': the handle of the figure you want to backup

% check if the file name exists
r = dir;
clear i
nametocheck = strcat(filename,filetype);
for i = 1 : length(r); 
    scanres(i) = strcmp(r(i).name,nametocheck); 
end

%if the file name exists, increment
i = 0;
if ~isempty(find(scanres==1))
    while ~isempty(find(scanres==1))
        i = i + 1;
        
        %increment the name
        filename = strcat(filename,'_0',num2str(i));
        nametocheck = strcat(filename,filetype);
        
        %check if the name exists already
        for j=1:length(r); scanres(j) = strcmp(r(j).name,nametocheck); end
    end
    errorms = [' "This filename exists already. The new filename is "',filename,'" '];
    disp(errorms)
end
if strcmp(filetype,'.fig')
    saveas(task,filename,'fig'); %name is a string
elseif strcmp(filetype,'.mat')
    save(filename,inputname(1)); %name is a string
end


