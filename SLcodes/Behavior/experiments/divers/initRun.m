
%initRun
%Date: 140418
%Author: Steeve laquitaine
%usage:
%[task]=initRun('steeve_fMRIexp01_metho_Pstd020_mean225_coh0061000_dir36_t103_029perCoh_140418',225,8.748,[.06 1],[103 29],5:10:355,[1 .2 0])

%requires my code library
%addpath(genpath('Users/steeve/DropBox/Library/Library_code'))


%Notes:
%a good example for our purpose: to choose a good std, we have to trade
%between:
%A)"setting a very narrow distribution and get prior's maximum effect, but
%keeping the full range of sample directions (lower probability directions may
%not occur at all").

%B)"setting the distribution width so that we always get at least 1 occurrence
%of the lower probability directions hence maintaining the full range of sample
%directions. But 1 occurence is not sufficient for statistics over the
%lower probability directions" --> no conclusions are possible for lower
%proba.directions

%C) "setting the distribution width so that we always get at least 5
%occurences of the lower proba.directions". The effect of the prior is reduced but
%we maintain a same full range of sample directions over conditions
%(cannot counfound). Moreover non parametric statistics becomes possible for lower probability
%directions.

%--> I chose A in this set of experiments to first check if we have any
%effect of the prior at all.

%plot the arrows on a unit circle
%modified from http://stackoverflow.com/questions/1803043/how-do-i-display-an-arrow-positioned-at-a-specific-angle-in-matlab


%Updates
%-------
%[121201]
%Current prior generation method is suceptible to noise due to random sampling
%e.g., distribution not peak at the mean or shape asymetry around the
%mean. I changed method. Now the number of repetition of each sample
%is approximated from a continuous distribution. This method ensures a
%perfectly symetric gaussian distribution peaking at the mean that is not
%susceptible to noise in random sampling.


%Inputs
%-------
%[121130]
%usage:
%You have to input the following parameters
%%Weak prior
%[task]=initPriors('steeve_exp07_metho_Pstd1000_mean225_card270_coh008_024_1_dir15_t60perCoh_121130',225,1000,60,155:10:295)
%note for a uniform distribution use:
%[task]=initPriors('steeve_exp07_metho_Pstdinf_mean225_card270_coh008_024_1_dir15_t60perCoh_121130',225,inf,60,155:10:295)
%%or Strong prior
%[task]=initPriors('steeve_exp07_metho_Pstd0020_mean225_card270_coh008_024_1_dir15_t60perCoh_121130',225,20,60,155:10:295)

%[121226]
%%Weak prior (uniform)
%[task]=initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,inf,[.12 .35 1],[105 75 30],155:10:295)
%%Strong prior (gaussian)
%[task]=initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,20,[.12 .35 1],[105 75 30],155:10:295)


%[130202]
%%No prior (uniform, all circle)
%[task]=initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,inf,[.12 .35 1],[105 75 30],0:10:350)

%%weak prior (slight bias)
%[task]=initRun('steeve_exp08_metho_Pstd100_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,100,[.12 .35 1],[105 75 30],155:10:295)

%%Intermediate weak prior (higher bias)
%[task]=initRun('steeve_exp08_metho_Pstd050_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,50,[.12 .35 1],[105 75 30],155:10:295)

%%clear prior (clear bias)
%[task]=initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,20,[.12 .35 1],[105 75 30],155:10:295)

%%Strong prior (strong bias)
%[task]=initRun('steeve_exp08_metho_Pstd14_mean225_card180270_coh012_035_1_dir15_t105_75_30perCoh_121206',225,14,[.12 .35 1],[105 75 30],155:10:295)


%[130204]
%%Strong prior
%[task]=initRun('steeve_exp08_metho_Pstd14_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, 15,[.08 .16 32],[105 75 30],5:10:355)

%%Less strong prior
%[task]=initRun('steeve_exp08_metho_Pstd20_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, 20,[.08 .16 32],[105 75 30],5:10:355)

%%Weaker prior
%[task]=initRun('steeve_exp08_metho_Pstd40_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, 40,[.08 .16 32],[105 75 30],5:10:355)

%%No prior
%[task]=initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh008_016_032_dir36_t144_108_36perCoh_130205',225,367,[.08 .16 .32],[144 108 36],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstdinf_mean225_card180270_coh008_016_032_dir36_t108_072_36perCoh_130205',225,367,[.08 .16 .32],[108  72 36],5:10:355)



%[130205]
%%Strong prior
%[task]=initRun('steeve_exp08_metho_Pstde2_50_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, exp(1.25*2),[.08 .16 .32],[108 72 36],5:10:355)

%%Less strong prior
%[task]=initRun('steeve_exp08_metho_Pstde3_75_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, exp(1.25*3),[.08 .16 .32],[108 72 36],5:10:355)

%%Weaker prior
%[task]=initRun('steeve_exp08_metho_Pstde5_00_mean225_card180270_coh008_016_032_dir36_t105_75_30perCoh_130205'  ,225, exp(1.25*4),[.08 .16 .32],[108 72 36],5:10:355)

%%No prior
%[task]=initRun('steeve_exp08_metho_Pstde6_25_mean225_card180270_coh008_016_032_dir36_t144_108_36perCoh_130205',225, exp(1.25*5),[.08 .16 .32],[144 108 36],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstde6_25_mean225_card180270_coh008_016_032_dir36_t108_072_36perCoh_130205',225, exp(1.25*5),[.08 .16 .32],[108  72 36],5:10:355)


%[130210]
%%Strong prior (even lower coherence)
%[task]=initRun('steeve_exp08_metho_Pstd012_mean225_card180270_coh004_012_028_dir36_t105_75_30perCoh_130209'  ,225, exp(1.25*2),[.04 .12 .28],[108 72 36],5:10:355)
%%Less strong prior
%[task]=initRun('steeve_exp08_metho_Pstd043_mean225_card180270_coh004_012_028_dir36_t105_75_30perCoh_130209'  ,225, exp(1.25*3),[.04 .12 .28],[108 72 36],5:10:355)
%%Weaker prior
%[task]=initRun('steeve_exp08_metho_Pstd149_mean225_card180270_coh004_012_028_dir36_t105_75_30perCoh_130209'  ,225, exp(1.25*4),[.04 .12 .28],[108 72 36],5:10:355)
%%No prior
%[task]=initRun('steeve_exp08_metho_Pstd518_mean225_card180270_coh004_012_028_dir36_t144_108_36perCoh_130209',225, exp(1.25*5),[.04 .12 .28],[144 108 36],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstd518_mean225_card180270_coh004_012_028_dir36_t108_072_36perCoh_130209',225, exp(1.25*5),[.04 .12 .28],[108  72 36],5:10:355)


%[130215]
%New priors (same coherence as previous)
%[task]=initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh004_012_028_dir36_t107_075_033perCoh_130215'  ,225, exp(1.25*2.4),[.04 .12 .28],[107 75 33],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstd020_mean225_card180270_coh008_016_032_dir36_t107_075_033perCoh_130215'  ,225, exp(1.25*2.4),[.08 .16 .32],[107 75 33],5:10:355)

%[task]=initRun('steeve_exp08_metho_Pstd016_mean225_card180270_coh004_012_028_dir36_t107_075_032perCoh_130215'  ,225, 16,[.04 .12 .28],[107 75 32],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstd016_mean225_card180270_coh008_016_032_dir36_t107_075_032perCoh_130215'  ,225, 16,[.08 .16 .32],[107 75 32],5:10:355)

%[task]=initRun('steeve_exp08_metho_Pstd018_mean225_card180270_coh004_012_028_dir36_t107_075_032perCoh_130215'  ,225, 18,[.04 .12 .28],[108 73 33],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstd018_mean225_card180270_coh008_016_032_dir36_t107_075_032perCoh_130215'  ,225, 18,[.08 .16 .32],[108 73 33],5:10:355)

%[task]=initRun('steeve_exp08_metho_Pstd030_mean225_card180270_coh004_012_028_dir36_t107_075_032perCoh_130215'  ,225, 30,[.04 .12 .28],[108 73 33],5:10:355)
%[task]=initRun('steeve_exp08_metho_Pstd030_mean225_card180270_coh008_016_032_dir36_t107_075_032perCoh_130215'  ,225, 30,[.08 .16 .32],[108 73 33],5:10:355)


%[130216]
%[task]=initRun('steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_073_033perCoh_130216', 225, 10,[.06 .12 .24],[107 73 33],5:10:355)
%[task]=initRun('steeve_exp12_metho_Pstd020_mean225_coh006012024_dir36_t107_075_033perCoh_130216', 225, 20,[.06 .12 .24],[107 75 33],5:10:355)
%[task]=initRun('steeve_exp12_metho_Pstd040_mean225_coh006012024_dir36_t107_075_033perCoh_130216', 225, 40,[.06 .12 .24],[101 75 31],5:10:355)
%[task]=initRun('steeve_exp12_metho_Pstd080_mean225_coh006012024_dir36_t107_075_033perCoh_130216', 225, 80,[.06 .12 .24],[106 75 34],5:10:355)

%[130217]
%[task]=initRun('steeve_exp12_metho_Pstd080_mean225_coh008_dir36_t107_075_033perCoh_130217', 225, 80,[.08],[106],5:10:355)

%[131224]
%0.74848 (std=80 degrees)
%2.7714 (std=40 degrees)
%8.748 (std=20 degrees)
%33.336 (std=10 degrees)
%[task]=initRun('steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_074_033perCoh_131224',225,33.336,[.06 .12 .24],[107 74 40],5:10:355,[.75 .75 0])
%[task]=initRun('steeve_exp12_metho_Pstd020_mean225_coh006012024_dir36_t103_074_029perCoh_131224',225,8.748,[.06 .12 .24],[103 74 29],5:10:355,[1 .2 0])
%[task]=initRun('steeve_exp12_metho_Pstd040_mean225_coh006012024_dir36_t100_075_032perCoh_131224',225,2.7714,[.06 .12 .24],[100 75 32],5:10:355,[1 .6 0])
%[task]=initRun('steeve_exp12_metho_Pstd080_mean225_coh006012024_dir36_t100_075_034perCoh_131224',225,0.74848,[.06 .12 .24],[110 74 40],5:10:355,[0.5 0 0])

%[140418]
%8.748 (std=20 degrees)
%[task]=initRun('steeve_fMRIexp01_metho_Pstd020_mean225_coh0061000_dir36_t103_029perCoh_140418',225,8.748,[.06 1],[103 29],5:10:355,[1 .2 0])


%a run
function task=initRun(Prname,Prmean,Prstd,coh,Prnumtrials_perCoh,PrSp,colorH)

%backup .m file in .txt file
filename=matlab.desktop.editor.getActiveFilename;
SLBackup(filename)

%Initialize figure
fig.hdle=figure('color','w');
fig.name=Prname;

%%%Prnumtrials_perCoh: is a vector indicating the number of trials of
%different conditions (e.g., coherences).
Prior.parameter.dir.trialnum=[];
%Prior.parameter.dir.series=[];
task.parameter.dir.series=[];
Priortmp={nan};
Prior.parameter.dir.count=[];
task.parameter.dir.coh=[];
for i=1:numel(Prnumtrials_perCoh)
    
    %Combine conditions of Priors (e.g., factor 1) and coherences (e.g., factor 2).
    [Priortmp{i}]=initPriors(Prname,Prmean,Prstd,Prnumtrials_perCoh(i),PrSp);
    
    %Initialize task per coherence
    task.parameter.dir.inputTrialnumperCoh(i)=Priortmp{i}.parameter.dir.trialnum;
    task.parameter.dir.TrueTrialnumperCoh(i)=numel(Priortmp{i}.parameter.dir.series);
    task.parameter.dir.trueStdperCoh(i)=Priortmp{i}.parameter.dir.trueStd;
    task.parameter.dir.countperCoh(i,:)=Priortmp{i}.parameter.dir.count;
    task.parameter.dir.seriesperCoh{i}=Priortmp{i}.parameter.dir.series;
    task.parameter.dir.pperCoh(i,:)=Priortmp{i}.parameter.dir.p;
    
    %Pool across coherences.
    task.parameter.dir.series=[task.parameter.dir.series Priortmp{i}.parameter.dir.series];
    repeat=repmat(coh(i),1,task.parameter.dir.TrueTrialnumperCoh(i));
    task.parameter.dir.coh=[task.parameter.dir.coh repeat];
end

%Pool across coherences.
task.parameter.dir.mean=Priortmp{1}.parameter.dir.mean;
task.parameter.dir.std=Priortmp{1}.parameter.dir.std;
task.parameter.dir.sample.degree=Priortmp{1}.parameter.dir.sample.degree;
task.parameter.dir.sampsiz=Priortmp{1}.parameter.dir.sampsiz;
task.parameter.dir.trialnum=sum(task.parameter.dir.TrueTrialnumperCoh);
task.parameter.dir.count=sum(task.parameter.dir.countperCoh);
task.parameter.dir.cohSample=coh;

%Draw circular histograms of motion directions for each coherence
for j=1:numel(coh)
    subplot(3,2,2*j-1);
    circHist(task.parameter.dir.sample.degree,task.parameter.dir.countperCoh(j,:),colorH,5)
end
fprintf('\n %12s \n','Drawing of circular histograms of motion directions for each coherence...done')

%Draw linear histograms of motion directions for each coherence
for j=1:numel(coh)
    subplot(3,2,j*2);
    
    %set plot's name
    Plotname=['coh ',num2str(coh(j))];
    
    %draw
    drawhistdir(task.parameter.dir.sample.degree,task.parameter.dir.countperCoh(j,:),fig,colorH);
    text(0.9*max(task.parameter.dir.sample.degree),.90*max(task.parameter.dir.countperCoh(j,:)),Plotname,'fontsize',14)
    xlabel({'Displayed directions','(in degree)'})
end
fprintf('\n %12s \n','Drawing of linear histograms of motion directions for each coherence...done')

%Backup figure and parameters
autobackup(fig.hdle,fig.name,'.fig')
autobackup(task,fig.name,'.mat')

%prior
function task=initPriors(Prname,Prmean,Prstd,Prnumtrials,PrSp)
%fig.name=Prname;
%mean, std of the distribution (stay the same).
task.parameter.dir.mean=Prmean;
task.parameter.dir.std=Prstd;
%nb of trials for each coherence in a sub-block
task.parameter.dir.trialnum=Prnumtrials;%65;
%you can choose to input or not a set of samples (if not, a default sample is chosen based on the sample size)
%task.parameter.dir.sample.degree=[90 100 110 120 130 140 150 160 170 180]; %set of sample directions (the artifacted directions in experiment 1 and 2)
task.parameter.dir.sample.degree=PrSp; %extend to include 90 degrees.
%task.parameter.dir.sample.degree=[40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220]; %extend to include cardinal directions
%nb of sample directions
task.parameter.dir.sampsiz=numel(task.parameter.dir.sample.degree);
%nb of occurence of each displayed directions
task.parameter.dir.count=[];

%if no series of directions has been input
if ~isfield(task.parameter.dir,{'series'})
    
    %if a set of sample direction has not been input either
    if ~isfield(task.parameter.dir,{'sample'})
        %create a default set of sample directions and a default series of
        %directions to be displayed on screen
        [task]=initSampleDir(task); %create a set of sample directions
        
        %Generate a discrete gaussian prior by random sampling (noisy).
        %[task]=initPriorDis_disc(task);
        %Generate a perfectly symetric gaussian prior
        [task]=initPriorDis_cont(task); %create a series of n trials directions
        disp(' "new directions have been initialized and drawn based on default parameters" ')
    else
        %if a set of displayed directions has been input.
        %Update the sample size
        task.parameter.dir.sampsiz=numel(task.parameter.dir.sample.degree);
        
        %Generate a discrete gaussian prior by random sampling (noisy).
        %[task]=initPriorDis_disc(task);
        %Generate a perfectly symetric gaussian prior
        [task]=initPriorDis_cont(task);
        disp(' "a series of n trials directions have been produced based on your input set of sample directions" ')
    end
    %Use input parameters.
else
    disp('the loaded directions statistics have been drawn')
end

%samples
function task=initSampleDir(task)

%initialize sample directions
for i=1:task.parameter.dir.sampsiz
    
    %calculate the directions(%Angle of arrow in degree, from x-axis)
    task.parameter.dir.sample.degree(i)=360/task.parameter.dir.sampsiz*i;
    
    %Rotate directions away from the cardinals.
    %    rotation=15; %in degree
    %or not
    rotation=0; %in degree
    task.parameter.dir.sample.degree(i)=task.parameter.dir.sample.degree(i) + rotation;
    
    %adjust the radian values according to the quadrant of the trigo.circle
    if task.parameter.dir.sample.degree(i)>360
        task.parameter.dir.sample.degree(i)=task.parameter.dir.sample.degree(i)-360;
    end
end
%sort the directions (adjusting the radian values creates an unsorted series of directions)
task.parameter.dir.sample.degree=sort(task.parameter.dir.sample.degree);

%continuous prior
function task=initPriorDis_cont(task)
%We make the probability densities we need based on von mises or gaussian
%laws. Then we calculate the number of occurrence of each direction
%(pdf*numTrials). Occurences are not real numbers because pdf are not but
%we need real numbers. We could randomly sample the density but
%for small sample size as it is the case in our experiment, this method
%creates small left/right asymetries relative to the mean. Those asymetries
%could create confounding biases in estimation that are not accounted for
%by the true continuous prior density.
%We avoid this asymetry by rounding the number of occurrence to real values
%and then by replicating direction trials according to the resulting number
%of occurrence. The incovenient is that the parameters of the prior std
%change but only slightly.
task.parameter.dir.series=[];
task.parameter.dir.count=[];
if task.parameter.dir.std~=inf
    
    %prior
    %PriorShape=gauss_distribution(task.parameter.dir.sample.degree,task.parameter.dir.mean,task.parameter.dir.std);
    PriorShape=vmPdfs(task.parameter.dir.sample.degree,task.parameter.dir.mean,task.parameter.dir.std,'norm');
    
    %create trials'count. We repeat samples according to the density.
    f2=PriorShape*task.parameter.dir.trialnum ;
    for i=1:numel(task.parameter.dir.sample.degree)
        f2repet=[];
        f2repet=repmat(task.parameter.dir.sample.degree(i),round(f2(i)),1);
        task.parameter.dir.series=[task.parameter.dir.series; f2repet];
    end
    task.parameter.dir.series=task.parameter.dir.series';
end

%count occurrence of displayed directions
for i=1:task.parameter.dir.sampsiz
    task.parameter.dir.count(i)=numel(find(task.parameter.dir.series==task.parameter.dir.sample.degree(i)));
end

%or generate a uniform prior
if task.parameter.dir.std==inf
    %check if "task.parameter.dir.trialnum" is a multiple of
    %"task.parameter.dir.sampsiz".
    if rem(task.parameter.dir.trialnum,task.parameter.dir.sampsiz)==0
        disp(['--- A UNIFORM distribution is being drawn ----'])
        task.parameter.dir.count=repmat(task.parameter.dir.trialnum/task.parameter.dir.sampsiz,...
            1,task.parameter.dir.sampsiz);
        task.parameter.dir.series=repmat(task.parameter.dir.sample.degree,task.parameter.dir.count(1),1);
        task.parameter.dir.series=task.parameter.dir.series(:);
        task.parameter.dir.series=task.parameter.dir.series';
    else
        disp(['"task.parameter.dir.trialnum" must be a multiple of "task.parameter.dir.sampsiz"'])
        return
    end
end

%calculate prior's actual std
task.parameter.dir.p=task.parameter.dir.count/sum(task.parameter.dir.count);
dis=task.parameter.dir.sample.degree-task.parameter.dir.mean;
task.parameter.dir.trueStd=sqrt(sum(task.parameter.dir.p.*(dis.^2)));
task.parameter.dir.trueMean=sum(task.parameter.dir.p.*task.parameter.dir.sample.degree);

%histogram
function axh=drawhistdir(sample,count,fig,colorH)
%graphics
axh.name=strcat(fig.name,'_his');
myBarWidth=1;

%draw
hold all

%histogram
bar_xend=2*sample(end)-sample(end-1);
bar([sample bar_xend],[count count(1)],...
    'FaceColor',colorH,...
    'BarWidth',myBarWidth,...
    'EdgeColor','none')

%lines
% plot([sample bar_xend],[count count(1)],...
%     'color',colorH,...
%     'lineWidth',2,...
%     'linesmoothing','on');
plot([225 225],[0 max(count)],'w:')

%Set x-axis
xunit=1:11:numel(sample);
set(gca,'fontsize',14,...
    'xtick',sample(xunit),...
    'xticklabel',sample(xunit))
ylabel('Count')
xlim([min(sample)-10 bar_xend+10])
ylim([0 max(count)]);
box off
%gaussian distribution
function functionf=gauss_distribution(x,mu,s)
%x:predicted variable
%mu:mean
%s:standard deviation
%height: it is the height of the gaussian peak
%p1=-.5 * ((x - mu)/s) .^ 2;

%before 131224....
%set the shaping parameters of the gaussian
%height=.5;
%height=.5;
%p1=-height*((x-mu)/s).^2;
%p2=(s*sqrt(2*pi));
%f=exp(p1)./p2;

%from 131224
f=(1/(s*sqrt(2*pi)))*exp(-((x-mu).^2)/(2*s.^2));
%based on http://en.wikipedia.org/wiki/Normal_distribution

%Convert from arbitrary unit to probabilities.
functionf=f/sum(f);
