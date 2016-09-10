

%SLinitDiscreteOrientationDis.m
%
%author: steeve laquitaine
%  date: 160111
%
% usage:
%
%
%         task = SLinitDiscreteOrientationDis('paramsPrior1',225,33.336,[1 0.12 0.1],[107 74 40],5:10:355,[.75 .75 0])
%         task = SLinitDiscreteOrientationDis('paramsPrior2',225,8.748,[1 0.12 0.1],[103 74 29],5:10:355,[1 .2 0])
%         task = SLinitDiscreteOrientationDis('paramsPrior3',225,2.7714,[1 0.12 0.1],[100 75 32],5:10:355,[1 .6 0])
%         task = SLinitDiscreteOrientationDis('paramsPrior4',225,0.74848,[1 0.12 0.1],[110 74 40],5:10:355,[0.5 0 0])
%
%Description:
%
% - at least 5 samples min of each direction (full range of sample dir).
%
% - symmetric distribution by discretization of a continuous dist. (instead of sampling
%   which produce assymetry)
%
% - 4 prior concentration params
%         0.74848 (std=80 degrees)
%         2.7714 (std=40 degrees)
%         8.748 (std=20 degrees)
%         33.336 (std=10 degrees)
%
%a run
function task = SLinitDiscreteOrientationDis(Prname,Prmean,Prstd,con,Prnumtrials_perCon,PrSp,colorH)

%Initialize figure
fig.hdle = figure('color','w');
fig.name = Prname;

%%%Prnumtrials_perCon: is a vector indicating the number of trials of
%different conditions (e.g., coherence or contrast).
Prior.parameter.loc.trialnum = [];
task.parameter.loc.series = [];
Priortmp  = {nan};
Prior.parameter.loc.count = [];
task.parameter.loc.con = [];
for i = 1 : numel(Prnumtrials_perCon)
    
    %Combine conditions of Priors (e.g., factor 1) and coherences (e.g., factor 2).
    [Priortmp{i}] = initPriors(Prname,Prmean,Prstd,Prnumtrials_perCon(i),PrSp);
    
    %Initialize task per coherence
    task.parameter.loc.inputTrialnumperCon(i) = Priortmp{i}.parameter.loc.trialnum;
    task.parameter.loc.TrueTrialnumperCon(i)  = numel(Priortmp{i}.parameter.loc.series);
    task.parameter.loc.trueStdperCon(i)       = Priortmp{i}.parameter.loc.trueStd;
    task.parameter.loc.countperCon(i,:)       = Priortmp{i}.parameter.loc.count;
    task.parameter.loc.seriesperCon{i}        = Priortmp{i}.parameter.loc.series;
    task.parameter.loc.pperCon(i,:)           = Priortmp{i}.parameter.loc.p;
    
    %Pool over coherences.
    task.parameter.loc.series = [task.parameter.loc.series Priortmp{i}.parameter.loc.series];
    repeat = repmat(con(i),1,task.parameter.loc.TrueTrialnumperCon(i));
    task.parameter.loc.con = [task.parameter.loc.con repeat];
end

%Pool over coherences.
task.parameter.loc.mean = Priortmp{1}.parameter.loc.mean;
task.parameter.loc.std = Priortmp{1}.parameter.loc.std;
task.parameter.loc.sample.degree = Priortmp{1}.parameter.loc.sample.degree;
task.parameter.loc.sampsiz = Priortmp{1}.parameter.loc.sampsiz;
task.parameter.loc.trialnum = sum(task.parameter.loc.TrueTrialnumperCon);
task.parameter.loc.count = sum(task.parameter.loc.countperCon);
task.parameter.loc.conSample = con;

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
function task = initPriors(Prname,Prmean,Prstd,Prnumtrials,PrSp)
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
        task = initPriorDis_cont(task);
        disp(' "a series of n trials directions have been produced based on your input set of sample directions" ')
    end
    %Use input parameters.
else
    disp('the loaded directions statistics have been drawn')
end

%samples
function task = initSampleDir(task)

%initialize sample directions
for i = 1 : task.parameter.loc.sampsiz
    
    %calculate the directions(%Angle of arrow in degree, from x-axis)
    task.parameter.loc.sample.degree(i) = 360/task.parameter.loc.sampsiz*i;
    
    %Rotate directions away from the cardinals.
    %    rotation=15; %in degree
    %or not
    rotation=0; %in degree
    task.parameter.loc.sample.degree(i) = task.parameter.loc.sample.degree(i) + rotation;
    
    %adjust the radian values according to the quadrant of the trigo.circle
    if task.parameter.loc.sample.degree(i)>360
        task.parameter.loc.sample.degree(i) = task.parameter.loc.sample.degree(i)-360;
    end
end
%sort the directions (adjusting the radian values creates an unsorted series of directions)
task.parameter.loc.sample.degree = sort(task.parameter.loc.sample.degree);

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
    PriorShape = vmPdfs(task.parameter.loc.sample.degree,task.parameter.loc.mean,task.parameter.loc.std,'norm');
    
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
function axh = drawhistdir(sample,count,fig,colorH)
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
%circular histogram
function SLcircHist(angles,lengths,colorH,baseRadius,iftext,varargin)
hold all

if ~exist('iftext','var')
    iftext=0;
end

%check that angles and lengths are column vectors
if size(angles,1)<size(angles,2)
    angles=angles';
end
if size(lengths,1)<size(lengths,2)
    lengths=lengths';
end

%don't show 0-length arrows
angles(lengths==0)=[];
lengths(lengths==0)=[];

%draw base circle

%calculate x0s & y0s
%We set the radius of the circle that supports the vectors in percent of
%the maximum vector length and draw the base circle.
%scale only if askedso.
if strcmp(varargin,'scaledRadius')
    baseRadius=baseRadius*max(lengths);
end

cAx=SLcircle(0,0,baseRadius,[0.75 0.75 0.75]);
axis square
axis off
c0=SLpolar2cartesian(angles,baseRadius);
x0=c0(:,1);
y0=c0(:,2);

%calculate xEnd & yEnd
c2=SLpolar2cartesian(angles,lengths);
cEnd=c0+c2;
xEnd=cEnd(:,1);
yEnd=cEnd(:,2);

%draw clean polar & bars
p=polar(0,max(lengths)+baseRadius);
h2=findall(gca,'type','patch');
delete(h2)
h=findall(gca,'type','line');
delete(intersect(h(2:end),p));
t=findall(gca,'type','text');
delete(t);
arrow([x0 y0],[xEnd yEnd],'Length',0,'BaseAngle',[],'Width',3,'facecolor',...
    colorH,'edgecolor',colorH)
rlim=max(lengths)+baseRadius;
axis([-1 1 -1 1]*rlim);

%text
if strcmp(iftext,'text')==1
    %lengths
    for i=1:numel(angles)
        if lengths(i)~=0
            text(1.12*xEnd(i),1.12*yEnd(i),num2str(lengths(i)),'fontsize',10,'HorizontalAlignment','center','verticalAlignment','middle')
        end
    end
end
%von Mises densities
function mPdfs = vmPdfs(x,u,k,type)


%check that x is a row vector
if size(x,1)>size(x,2)
    x=x';
end

%check that u and k are col vectors
if size(u,1)<size(u,2)
    u=u';
end
%check that u is a col vector
if size(k,1)<size(k,2)
    k=k';
end

%radians
xrad = SLde2r(x,1); xrad=xrad';
urad = SLde2r(u,1); urad=urad';
k=k';
%When von mises with different mean u1,u2,u3 but with same k are input
%We can get von Mises with mean u2,u3,etc...simply by rotating the von
%mises with mean u1 by u2-u1, u3-u1 etc...
%When we don't do that we get slightly different von mises with different
%peakvalue due to numerical instability caused by cosine and exponential
%functions.
%case all k are same
if sum(k - k(1))==0
    %if mean u is not one of x
    if isempty(intersect(x,u)==0)
        fprintf('\n %s \n','(vmPdfs) WARNING : Mean has to be one of x values. Change mean.....')
        if isnan(u)
            sprintf('(vmPdfs) WARNING : Von Mises Mean is NaN....')
        end
        dbstack
        keyboard
    else
        
        %-------------------------
        %case k tends toward + inf
        %-------------------------
        %delta function
        if k(1) > 1e300
            
            mPdfs = zeros(length(xrad),1);
            mPdfs(xrad == urad(1)) = 1;
            %make other densities by circshifting the first
            for i = 2 : numel(u)
                rotation = find(x==u(i)) - find(x==u(1));
                mPdfs(:,i) = circshift(mPdfs(:,1),[rotation,0]);
            end
            
        else
            
            %---------
            %k non inf
            %---------
            w = 1;
            k = k(1);
            
            %should work with parallel computing (parfor loop). ".*" seems to
            %work badly with parfor. But I think k must be unique.
            mPdfs = exp(k*cos(w*(xrad-urad(1)))-k)./(2*pi.*besseli(0,k,1));
            %make other densities by circshifting the first
            for i = 2 : numel(u)
                rotation = find(x==u(i))-find(x==u(1));
                mPdfs(:,i) = circshift(mPdfs(:,1),[rotation,0]);
            end
            
            %scale to probabilities.
            if strcmp(type,'norm')==1
                %mPdfs=mPdfs./sum(mPdfs(:,1));
                mPdfs=mPdfs/sum(mPdfs(:,1));
            end
        end
    end
end

%When von mises with different mean and k are input, calculate densities
%separately. Make sure each mean u maps with a k.
if sum(k - k(1))~=0
    
    %case k not too high
    k     = k(ones(numel(xrad),1),:);
    x2rad = xrad(:,ones(numel(urad),1));
    u2rad = urad(ones(numel(xrad),1),:);
    k2    = k(ones(numel(xrad),1),:);
    w     = 1;
    mPdfs = exp(k2.*cos(w*(x2rad-u2rad))-k2)./(2*pi.*besseli(0,k2,1));
    
    %case k tends toward inf, delta density
    kinfCols      = find(k2(1,:)>300);
    u2radkinfCols = u2rad(:,kinfCols);
    deltas        = zeros(numel(xrad),length(kinfCols));
    deltas(x2rad(:,kinfCols)-u2radkinfCols==0) = 1;
    mPdfs(:,kinfCols) = deltas;
    
    %for parallel processing (bu 30 times slower)
    %     for i=1:size(k2,1)
    %         for j=1:size(k2,2)
    %             mPdfs(i,j)=exp(k2(i,j)*cos(w*(x2rad(i,j)-u2rad(i,j)))-k2(i,j))/(2*pi*besseli(0,k2(i,j),1));
    %         end
    %     end
    
    %scale to pdfs.
    if strcmp(type,'norm')==1
        Z_=sum(mPdfs);
        Z=Z_(ones(numel(x),1),:);
        mPdfs=mPdfs./Z;
    else
    end
end




