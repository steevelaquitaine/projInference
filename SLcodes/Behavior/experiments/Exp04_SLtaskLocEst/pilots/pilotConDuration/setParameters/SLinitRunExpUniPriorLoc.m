

%[121126] - Steeve
%
% usage
%
%INPUT
    %mean
    %variance
    %trialnum

%OUTPUT
    %vector of n trials directions: task.parameter.dir.series

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

%priors
%0.74848 (std=80 degrees)
%2.7714 (std=40 degrees)

% task = SLinitRunExpUniPriorLoc('steeve_exp04_metho_Pstd040_mean225_con000090000200003_loc36_t100_075_032perCon_151107',225,2.7714,[0.0009 0.0018 0.0036],[100 75 32],5:10:355,[1 .6 0])
% task = SLinitRunExpUniPriorLoc('steeve_exp04_metho_Pstd080_mean225_con000010000200003_loc36_t100_075_034perCon_151107',225,0.74848,[0.0009 0.0018 0.0036],[110 74 40],5:10:355,[0.5 0 0])
% task = SLinitRunExpUniPriorLoc('steeve_exp04_metho_Pstd040_mean225_con00014_00027_00054_loc36_t100_075_032perCon_151107',225,2.7714,[0.0014 0.0027 0.0054],[100 75 32],5:10:355,[1 .6 0])
% task = SLinitRunExpUniPriorLoc('steeve_exp04_metho_Pstd080_mean225_con00014_00027_00054_loc36_t100_075_034perCon_151107',225,0.74848,[0.0014 0.0027 0.0054],[110 74 40],5:10:355,[0.5 0 0])


%a run
function task=SLinitRunExpUniPriorLoc(Prname,Prmean,Prstd,con,Prnumtrials_perCon,PrSp,colorH)
%Initialize figure
fig.hdle=figure('color','w');
fig.name=Prname; 

%%%Prnumtrials_perCon: is a vector indicating the number of trials of
%different conditions (e.g., conerences).
Prior.parameter.loc.trialnum=[];
% Prior.parameter.loc.series=[];
task.parameter.loc.series=[];
Priortmp={nan};
Prior.parameter.loc.count=[];
task.parameter.loc.con=[];
for i=1:numel(Prnumtrials_perCon)
    
    %Combine conditions of Priors (e.g., factor 1) and coherences (e.g., factor 2).
    [Priortmp{i}]=initPriors(Prname,Prmean,Prstd,Prnumtrials_perCon(i),PrSp);
    
    %Initialize task per coherence
    task.parameter.loc.inputTrialnumperCon(i)=Priortmp{i}.parameter.loc.trialnum;
    task.parameter.loc.TrueTrialnumperCon(i)=numel(Priortmp{i}.parameter.loc.series);
    task.parameter.loc.trueStdperCon(i)=Priortmp{i}.parameter.loc.trueStd;
    task.parameter.loc.countperCon(i,:)=Priortmp{i}.parameter.loc.count;
    task.parameter.loc.seriesperCon{i}=Priortmp{i}.parameter.loc.series;
    task.parameter.loc.pperCon(i,:)=Priortmp{i}.parameter.loc.p;
    
    %Pool across coherences.
    task.parameter.loc.series=[task.parameter.loc.series Priortmp{i}.parameter.loc.series];
    repeat=repmat(con(i),1,task.parameter.loc.TrueTrialnumperCon(i));
    task.parameter.loc.con=[task.parameter.loc.con repeat];
end

%Pool across coherences.
task.parameter.loc.mean=Priortmp{1}.parameter.loc.mean;
task.parameter.loc.std=Priortmp{1}.parameter.loc.std;
task.parameter.loc.sample.degree=Priortmp{1}.parameter.loc.sample.degree;
task.parameter.loc.sampsiz=Priortmp{1}.parameter.loc.sampsiz;
task.parameter.loc.trialnum=sum(task.parameter.loc.TrueTrialnumperCon);
task.parameter.loc.count=sum(task.parameter.loc.countperCon);
task.parameter.loc.conSample=con;

%Draw distributions as polar plots
subplot(3,2,1);
SLcircHist(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(1,:),colorH,5)
subplot(3,2,3);
SLcircHist(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(2,:),colorH,5)
subplot(3,2,5);
SLcircHist(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(3,:),colorH,5)

%Draw distributions as histograms (con)
subplot(3,2,2);
drawhistdir(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(1,:),fig,colorH);
text(0.9*max(task.parameter.loc.sample.degree),.90* ...
     max(task.parameter.loc.countperCon(1,:)),['con' num2str(con(1))],'fontsize',14)
subplot(3,2,4);
drawhistdir(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(2,:),fig,colorH);
text(0.9*max(task.parameter.loc.sample.degree),.90*max(task.parameter.loc.countperCon(2,:)),['con' num2str(con(2))],'fontsize',14)
subplot(3,2,6);
drawhistdir(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(3,:),fig,colorH);
text(0.9*max(task.parameter.loc.sample.degree),.90*max(task.parameter.loc.countperCon(3,:)),['con' num2str(con(3))],'fontsize',14)
xlabel({'Displayed directions','(in degree)'})

%Backup figure
autobackup(fig.hdle,fig.name,'.fig')

%Backup parameters
autobackup(task,fig.name,'.mat')

%prior
function task=initPriors(Prname,Prmean,Prstd,Prnumtrials,PrSp)
%fig.name=Prname; 
%mean, std of the distribution (stay the same).
task.parameter.loc.mean=Prmean;   
task.parameter.loc.std=Prstd;
%nb of trials for each coherence in a sub-block
task.parameter.loc.trialnum=Prnumtrials;%65; 
%you can choose to input or not a set of samples (if not, a default sample is chosen based on the sample size) 
%task.parameter.loc.sample.degree=[90 100 110 120 130 140 150 160 170 180]; %set of sample directions (the artifacted directions in experiment 1 and 2) 
task.parameter.loc.sample.degree=PrSp; %extend to include 90 degrees.
%task.parameter.loc.sample.degree=[40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220]; %extend to include cardinal directions
%nb of sample directions
task.parameter.loc.sampsiz=numel(task.parameter.loc.sample.degree);
%nb of occurence of each displayed directions
task.parameter.loc.count=[];

%if no series of directions has been input
if ~isfield(task.parameter.loc,{'series'})
    
    %if a set of sample direction has not been input either
    if ~isfield(task.parameter.loc,{'sample'})    
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
        task.parameter.loc.sampsiz=numel(task.parameter.loc.sample.degree);
        
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
for i=1:task.parameter.loc.sampsiz
    
    %calculate the directions(%Angle of arrow in degree, from x-axis)
    task.parameter.loc.sample.degree(i)=360/task.parameter.loc.sampsiz*i;         
    
    %Rotate directions away from the cardinals.
%    rotation=15; %in degree
    %or not
    rotation=0; %in degree
    task.parameter.loc.sample.degree(i)=task.parameter.loc.sample.degree(i) + rotation;

    %adjust the radian values according to the quadrant of the trigo.circle
    if task.parameter.loc.sample.degree(i)>360
       task.parameter.loc.sample.degree(i)=task.parameter.loc.sample.degree(i)-360;
    end
end
%sort the directions (adjusting the radian values creates an unsorted series of directions)
task.parameter.loc.sample.degree=sort(task.parameter.loc.sample.degree);

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
task.parameter.loc.series=[];
task.parameter.loc.count=[];
if task.parameter.loc.std~=inf
    
    %prior
    %PriorShape=gauss_distribution(task.parameter.loc.sample.degree,task.parameter.loc.mean,task.parameter.loc.std);
    PriorShape=vmPdfs(task.parameter.loc.sample.degree,task.parameter.loc.mean,task.parameter.loc.std,'norm');

    %create trials'count. We repeat samples according to the density.
    f2=PriorShape*task.parameter.loc.trialnum ;
    for i=1:numel(task.parameter.loc.sample.degree)
        f2repet=[];
        f2repet=repmat(task.parameter.loc.sample.degree(i),round(f2(i)),1);
        task.parameter.loc.series=[task.parameter.loc.series; f2repet];
    end
    task.parameter.loc.series=task.parameter.loc.series';
end

%count occurrence of displayed directions
for i=1:task.parameter.loc.sampsiz
    task.parameter.loc.count(i)=numel(find(task.parameter.loc.series==task.parameter.loc.sample.degree(i)));
end
    
%or generate a uniform prior
if task.parameter.loc.std==inf
    %check if "task.parameter.loc.trialnum" is a multiple of
    %"task.parameter.loc.sampsiz".
    if rem(task.parameter.loc.trialnum,task.parameter.loc.sampsiz)==0
        disp(['--- A UNIFORM distribution is being drawn ----'])
        task.parameter.loc.count=repmat(task.parameter.loc.trialnum/task.parameter.loc.sampsiz,...
            1,task.parameter.loc.sampsiz);
        task.parameter.loc.series=repmat(task.parameter.loc.sample.degree,task.parameter.loc.count(1),1);
        task.parameter.loc.series=task.parameter.loc.series(:);
        task.parameter.loc.series=task.parameter.loc.series';
    else
        disp(['"task.parameter.loc.trialnum" must be a multiple of "task.parameter.loc.sampsiz"'])
        return
    end
end

%calculate prior's actual std
task.parameter.loc.p=task.parameter.loc.count/sum(task.parameter.loc.count);
dis=task.parameter.loc.sample.degree-task.parameter.loc.mean;
task.parameter.loc.trueStd=sqrt(sum(task.parameter.loc.p.*(dis.^2)));
task.parameter.loc.trueMean=sum(task.parameter.loc.p.*task.parameter.loc.sample.degree);

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
%backup figures&parameters
function autobackup(task,filename,filetype)
%input:
%structure called 'fig' with 2 fields:
%- 'name': e.g., '120910_Steeve_exp02_resul_run8_var1000_mean195_con024_t250_sess2'
%- 'hdle': the handle of the figure you want to backup

%check if the file name exists
r=dir;
clear i
nametocheck=strcat(filename,filetype);
for i=1:length(r); scanres(i)=strcmp(r(i).name,nametocheck); end

%if the file name exists, increment
i=0;
if ~isempty(find(scanres==1))
    while ~isempty(find(scanres==1))
        i=i+1;
        %increment the name
        filename=strcat(filename,'_0',num2str(i));
        nametocheck=strcat(filename,filetype);
        %check if the name exists already
        for j=1:length(r); scanres(j)=strcmp(r(j).name,nametocheck); end
    end
    errorms=[' "This filename exists already. The new filename is "',filename,'" '];
    disp(errorms)
end
if strcmp(filetype,'.fig')
    saveas(task,filename,'fig'); %name is a string
elseif strcmp(filetype,'.mat')
    save(filename,inputname(1)); %name is a string
end


