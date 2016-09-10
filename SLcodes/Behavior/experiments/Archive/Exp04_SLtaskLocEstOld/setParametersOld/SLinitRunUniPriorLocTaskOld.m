
%SLinitRunUniPriorLocTask.m
%
%  author: steeve laquitaine
%    date: 2015/04/22
% purpose: initialize the task parameters of a scan run of the location
%          task. Location vector is obtained from randomly drawing from 
%           a von Mises distribution. Contrasts trials are randomly
%           interleaved.
%          The ouput "task.mat" is used as input to SLtaskLocEst.m
%
% usage:
%   
%   80 deg prior
%   task = SLinitRunUniPriorLocTask('steeve_metho_Pstd080_mean225_con001100_dir06_t160_50perCon_training_150501',225,0.74848,[.01 1],[160 50],45:60:355,[.5 0 0]);
%
%   40º prior
%   task = SLinitRunUniPriorLocTask('steeve_metho_Pstd040_mean225_con001100_dir06_t160_51perCon_training_150501',225,2.7714 ,[.01 1],[160 50],45:60:355,[1 0.6 0]);
%
% non optional inputs:
%--------------------
%
%                   Prname : output file name
%                            e.g., 'steeve_exp12_metho_Pstd010_mean225_con006012024_dir36_t107_074_033perCon_131224'
%                   Prmean : mean of the von Mises distribution 
%                            e.g.,225
%                    Prstd : von Mises k parameter of the von Mises distribution
%                            e.g., 33.336 (std of 10 deg)
%                      con : contrast parameter(s)
%                            e.g., [0.1 0.5 1]
%       Prnumtrials_perCon : number of trials with each contrast parameter
%                            e.g., [107 74 40]
%                     PrSp : samples of the von Mises distribution
%                            e.g., 5:10:355     
%                   colorH : color of the plot of the von Mises distribution
%                            e.g., [.75 .75 0]   
%note:
%       k=0.74848 (std=80 degrees)
%       k=2.7714 (std=40 degrees)
%       k=8.748 (std=20 degrees)
%       k=33.336 (std=10 degrees)
%

function task = SLinitRunUniPriorLocTask(Prname,Prmean,Prstd,con,Prnumtrials_perCon,PrSp,colorH)

%initialize task parameters
task.parameter.loc.numCon = size(con,2);
task.parameter.loc.strength = Prstd;

%figure
fig.hdle = figure('color','w');
fig.name = Prname; 

%%%Prnumtrials_perCon: is a vector indicating the number of trials of
%different conditions (e.g., contrasts).
Prior.parameter.loc.trialnum=[];
% Prior.parameter.loc.series=[];
task.parameter.loc.series=[];
Priortmp={nan};
Prior.parameter.loc.count=[];
task.parameter.loc.con=[];
for i=1:numel(Prnumtrials_perCon)
    
    %Combine conditions of Priors (e.g., factor 1) and contrasts (e.g., factor 2).
    Priortmp{i} = initPriors(Prname,Prmean,Prstd,Prnumtrials_perCon(i),PrSp);
    
    %init task parameters for each contrast
    task.parameter.loc.inputTrialnumperCon(i) = Priortmp{i}.parameter.loc.trialnum;
    task.parameter.loc.TrueTrialnumperCon(i) = numel(Priortmp{i}.parameter.loc.series);
    task.parameter.loc.trueStdperCon(i) = Priortmp{i}.parameter.loc.trueStd;
    task.parameter.loc.countperCon(i,:) = Priortmp{i}.parameter.loc.count;
    task.parameter.loc.seriesperCon{i} = Priortmp{i}.parameter.loc.series;
    task.parameter.loc.pperCon(i,:) = Priortmp{i}.parameter.loc.p;
    task.parameter.loc.contDistInCount(i,:) = Priortmp{i}.parameter.loc.contDistInCount;
    
    %Pool across contrasts.
    task.parameter.loc.series=[task.parameter.loc.series Priortmp{i}.parameter.loc.series];
    repeat=repmat(con(i),1,task.parameter.loc.TrueTrialnumperCon(i));
    task.parameter.loc.con=[task.parameter.loc.con repeat];
end

%Pool across contrasts.
task.parameter.loc.mean = Priortmp{1}.parameter.loc.mean;
task.parameter.loc.std = Priortmp{1}.parameter.loc.std;
task.parameter.loc.sample.degree = Priortmp{1}.parameter.loc.sample.degree;
task.parameter.loc.sampsiz = Priortmp{1}.parameter.loc.sampsiz;
task.parameter.loc.trialnum = sum(task.parameter.loc.TrueTrialnumperCon);
task.parameter.loc.count = sum(task.parameter.loc.countperCon);
task.parameter.loc.conSample = con;

%polar and hist distributions
%Draw polar distributions
myOdd = SLmakeOdd(0,size(con,2)*2);
for i = 1 : size(con,2)
    subplot(size(con,2),2,myOdd(i));
    SLcircHist(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(i,:),colorH,5)
end

%Draw discrete distributions (linear space)
myEven = SLmakeEven(0,size(con,2)*2);
for i = 1: size(con,2)
    subplot(size(con,2),2,myEven(i));
    
    %discrete distribution
    SLdrawhistdir(task.parameter.loc.sample.degree,task.parameter.loc.countperCon(i,:),fig,colorH);
    
    %overlap continuous von Mises
    plot(task.parameter.loc.sample.degree,task.parameter.loc.contDistInCount(i,:),'color','k');
    text(0.9*max(task.parameter.loc.sample.degree),.90*max(task.parameter.loc.countperCon(i,:)),['con ' num2str(con(i))],'fontsize',14)
    
    %info
    title(['Prior k : ' num2str(Prstd),', trial nb : ' num2str(task.parameter.loc.TrueTrialnumperCon(i))])
    lg=legend('Actual disc. dis','','Orig.cont.dis','location','NorthWest');
    set(lg,'fontsize',10)
    legend('boxoff')
    
    %graphics
    maxp = max([task.parameter.loc.countperCon(i,:) task.parameter.loc.contDistInCount(i,:)]);
    ylim([0 maxp])
end
xlabel({'Displayed locations (deg)'})

%clean up
% SLConventionUp

%Backup figure
autobackup(fig.hdle,fig.name,'.fig')

%Backup parameters
autobackup(task,fig.name,'.mat')

%prior
function task = initPriors(Prname,Prmean,Prstd,Prnumtrials,PrSp)

%mean, std of the distribution (stay the same).
task.parameter.loc.mean = Prmean;   
task.parameter.loc.std = Prstd;

%nb of trials for each contrast in a sub-block
task.parameter.loc.trialnum = Prnumtrials;%65; 

%you can choose to input or not a set of samples 
%(if not, a default sample is chosen based on the sample size) 
task.parameter.loc.sample.degree = PrSp;

%nb of sample directions
task.parameter.loc.sampsiz=numel(task.parameter.loc.sample.degree);

%nb of occurence of each displayed directions
task.parameter.loc.count = [];

%if no series of directions has been input
if ~isfield(task.parameter.loc,{'series'})
    
    %if no input sample directions
    if ~isfield(task.parameter.loc,{'sample'})    
        
        %create a default set of sample directions and a default series of
        %directions to be displayed on screen
        %create a set of sample directions
        task = initSampleDir(task);
        
        %Generate a discrete Gaussian prior by random sampling (noisy).
        %task = initPriorDis_disc(task);
        
        %Generate a perfectly symetric Gaussian prior
        %create a series of n trials directions
        task = initPriorDis_cont(task); 
        fprintf('\n (SLinitRunUniPriorLocTask) "I have initialized new series of parameters and plotted based on default parameters" \n')
    else
        %if displayed directions have been input update the sample size
        task.parameter.loc.sampsiz = numel(task.parameter.loc.sample.degree);
        
        %Generate a discrete Gaussian prior by random sampling (noisy).
        %task = initPriorDis_disc(task);
        %Generate a perfectly symetric Gaussian prior
        task = initPriorDis_cont(task);
        fprintf(' \n (SLinitRunUniPriorLocTask) "I have produced a series of n. directions based on your input sample directions" \n')
    end
    %Use input parameters.
else
    fprintf('\n (SLinitRunUniPriorLocTask) I have plotted the loaded directions statistics \n')
end
 
%samples
function task = initSampleDir(task)

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
function task = initPriorDis_cont(task)
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
    %PriorShape = gauss_distribution(task.parameter.loc.sample.degree,task.parameter.loc.mean,task.parameter.loc.std);
    PriorShape = vmPdfs(task.parameter.loc.sample.degree,task.parameter.loc.mean,task.parameter.loc.std,'norm');

    %create trials' count. We repeat samples according to the density.
    contDistInCount = PriorShape*task.parameter.loc.trialnum ;
    for i = 1 : numel(task.parameter.loc.sample.degree)
        contDistInCountrepet = [];
        contDistInCountrepet = repmat(task.parameter.loc.sample.degree(i),round(contDistInCount(i)),1);
        task.parameter.loc.series = [task.parameter.loc.series; contDistInCountrepet];
    end
    task.parameter.loc.series = task.parameter.loc.series';
end

%count occurrence of displayed directions
for i=1:task.parameter.loc.sampsiz
    task.parameter.loc.count(i) = numel(find(task.parameter.loc.series==task.parameter.loc.sample.degree(i)));
end
    
%or generate a uniform prior
if task.parameter.loc.std==inf
    %check if "task.parameter.loc.trialnum" is a multiple of
    %"task.parameter.loc.sampsiz".
    if rem(task.parameter.loc.trialnum,task.parameter.loc.sampsiz)==0
        fprintf(['--- A UNIFORM distribution is being drawn ----'])
        task.parameter.loc.count=repmat(task.parameter.loc.trialnum/task.parameter.loc.sampsiz,...
            1,task.parameter.loc.sampsiz);
        task.parameter.loc.series=repmat(task.parameter.loc.sample.degree,task.parameter.loc.count(1),1);
        task.parameter.loc.series=task.parameter.loc.series(:);
        task.parameter.loc.series=task.parameter.loc.series';
    else
        fprintf(['\n (SLinitRunUniPriorLocTask) "task.parameter.loc.trialnum" must be a multiple of "task.parameter.loc.sampsiz" \n'])
        return
    end
end

%calculate prior's actual std
task.parameter.loc.p = task.parameter.loc.count/sum(task.parameter.loc.count);
dis = task.parameter.loc.sample.degree-task.parameter.loc.mean;
task.parameter.loc.trueStd = sqrt(sum(task.parameter.loc.p.*(dis.^2)));
task.parameter.loc.trueMean = sum(task.parameter.loc.p.*task.parameter.loc.sample.degree);
task.parameter.loc.contDistInCount = contDistInCount;

%histogram
% function axh = drawhistdir(sample,count,fig,colorH)
% %graphics
% axh.name=strcat(fig.name,'_his'); 
% myBarWidth=1;
% 
% %draw
% hold all
% 
% %histogram
% bar_xend=2*sample(end)-sample(end-1);
% bar([sample bar_xend],[count count(1)],...
%     'FaceColor',colorH,...
%     'BarWidth',myBarWidth,...
%     'EdgeColor','none')
% 
% %lines
% % plot([sample bar_xend],[count count(1)],...
% %     'color',colorH,...
% %     'lineWidth',2,...
% %     'linesmoothing','on');
% plot([225 225],[0 max(count)],'w:')
% 
% %annotate
% for i = numel()
% text([sample bar_xend],[count count(1)],num2str(count(j)),'fontsize',10)
% 
% %Set x-axis
% xunit=1:11:numel(sample);
% set(gca,'fontsize',14,...
%     'xtick',sample(xunit),...
%     'xticklabel',sample(xunit))
% ylabel('Count')
% xlim([min(sample)-10 bar_xend+10])
% ylim([0 max(count)]);
% box off



%gaussian distribution
function functionf = gauss_distribution(x,mu,s)
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
r = dir;
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
    fprintf(['\n (SLinitRunUniPriorLocTask) This filename exists already. I created a new file named "',filename,'" \n'])
end
if strcmp(filetype,'.fig')
    saveas(task,filename,'fig'); %name is a string
elseif strcmp(filetype,'.mat')
    save(filename,inputname(1)); %name is a string
end


