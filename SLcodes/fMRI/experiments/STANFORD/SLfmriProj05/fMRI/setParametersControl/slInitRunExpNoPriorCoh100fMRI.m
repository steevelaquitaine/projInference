
%slInitRunExpNoPriorCoh100fMRI.m 
%[121126] - Steeve
%
% usage
%
% INPUT
%     mean
%     variance
%     trialnum
% 
% OUTPUT
%     vector of n trials directions: task.parameter.dir.series
% 
% 
% task = slInitRunExpNoPriorCoh100fMRI('fMRIParamsNoPriorCoh100',225,0,1,257,15:70:355,[0.5 0 0])

%a run
function task = slInitRunExpNoPriorCoh100fMRI(Prname,Prmean,Prstd,coh,Prnumtrials_perCoh,PrSp,colorH)

%Initialize figure
fig.hdle = figure('color','w');
fig.name = Prname; 

%%%Prnumtrials_perCoh: is a vector indicating the number of trials of
%different conditions (e.g., coherences).
Prior.parameter.dir.trialnum=[];
% Prior.parameter.dir.series=[];
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

%Draw distributions as polar plots
nCoh = length(coh);
pos = [1 : nCoh] + [0 : nCoh-1];
for i = 1 : nCoh
    subplot(nCoh,2,pos(i));
    SLcircHist(task.parameter.dir.sample.degree,task.parameter.dir.countperCoh(i,:),colorH,5)
end

%Draw distributions as histograms (coh)
pos = [1 : nCoh]*2;
for i = 1 : nCoh
    subplot(nCoh,2,pos(i));
    drawhistdir(task.parameter.dir.sample.degree,task.parameter.dir.countperCoh(i,:),fig,colorH);
    
    %overlap continuous von Mises
    text(0.9*max(task.parameter.dir.sample.degree),.90*max(task.parameter.dir.countperCoh(i,:)),['con ' num2str(coh(i))],'fontsize',14)
    
    %info
    for j = 1 : task.parameter.dir.sampsiz
        y = task.parameter.dir.countperCoh(i,j)/2;
        mytxt = num2str(task.parameter.dir.countperCoh(i,j));
        text(task.parameter.dir.sample.degree(j),y,mytxt)
    end
    title(['Prior k : ' num2str(Prstd),', trial nb : ' num2str(task.parameter.dir.TrueTrialnumperCoh(i))])
    lg = legend('Actual disc. dis','','Orig.cont.dis','location','NorthWest');
    set(lg,'fontsize',10)
    legend('boxoff')
     
    %graphics
    maxp = max(task.parameter.dir.countperCoh(i,:));
    ylim([0 maxp])
    set(gca,'xtick',task.parameter.dir.sample.degree(1:1:end),'xticklabel',task.parameter.dir.sample.degree(1:1:end))
    
    
end
xlabel({'Displayed directions','(in degree)'})

%Backup figure
autobackup(fig.hdle,fig.name,'.fig')

%Backup parameters
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
xunit=1 : 1 : numel(sample);
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
%- 'name': e.g., '120910_Steeve_exp02_resul_run8_var1000_mean195_coh024_t250_sess2'
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


