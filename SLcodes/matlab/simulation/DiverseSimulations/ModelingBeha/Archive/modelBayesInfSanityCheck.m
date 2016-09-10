

%....................
%Features of the code
%....................
%Data are linearized by adding or substracting 360 if needed. Thus
%modeling is possible.

%...........................
%The models used in the code
%...........................
%Currently we made the assumption that the noise in the sensory evidence (std(LLH))
%is represented linearly.
%(see  Hurlimann,F.,Kiper,D. C.,& Carandini,M. (2002).
%Testing the Bayesian model of perceived speed. Vision Res,42(19),2253?2257.).

%Different assumptions are made concerning the representation of the noise
%in the sensory evidence (std(LLH))

%#Assumption 1-model 1.
%std LLH=Sreal * Co^n/(Co^n + C50^n) with n=2 if the representation of
%coherence is linearly encoded.

%#Assumption 2-model 2 (simplest)
%std LLH=Sreal * C

%#Assumption 2-model 2 (simplest)
%std LLH=Sreal * C + So (So ~=0)


%.......
%Updates
%.......
%130415 - Normalization of the circular data estimated direction=
%distance (displayed,estimated) + displayed.
%This manipulation linearizes the data.

%130419 -  Fit the model to mean and variance of the data pooled together
%I cannot derive the equations for the posterior variance and mean such that I
%get only one free parameter.
%I use the task design (independent manipulation of the prior's variance
%and likelihood's variance) to fit the model with two free parameters:
%prior's variance and likelihood'variance.

%.....
%To do
%.....
%130511 - clean up and improve BI3-5....
%130516 - add descriptive statistics


%bug


%Models
function [Sp,pred,fitP,R2,udata,sdata,dataOut,predOut,Fbp,F,stdPa,fitPt,logl]=modelBayesInfSanityCheck(data,Fs,fig,FofLl2del,Ll2del)
pred=[];
fitP=[];
R2=[];
udata=[];
sdata=[];
dataOut=[];
predOut=[];
Fbp=[];
F=[];
Sp=[];
stdPa=[];
fitPt=[];

%databank
tic
fprintf('\n %12s \n','Now create data matrix...')
Makedatabank(data,FofLl2del,Ll2del);
toc

%factors
tic
fprintf('\n %12s \n','Now sorting the factors...')
[datafitm,datafitt,F]=initFs(Fs,fig);
toc

%.........................................................................
%model 1
% Fit data's mean with width of llh(std) an inverse function of coherence.
%initial parameters for k were chosen such that simulated uPo span
%different situations ranging from totally biased toward the likelihood or
%totally biased toward the prior.For visualization run the following.
% uLl=5:15:355;
% for i=1:numel(uLl)
%     k=0:40:119; uPr=225; c=6; uPo=[(1./(1+(k./c).^2)).*uLl(i)+uPr.*(1./(1+(c./k).^2))]';
%     hold on; plot(uPo);
% end
%min/max fitting range was chosen between 0 and +inf where the model's
%equation can be solved.
%You are free to set (or not,[]) 1 free parameters.
%.........................................................................
% fitP=[];%[0.1 0.2 0.3 0.4];
% [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa]=LSft(datafitm,fitP);
% drawPred(pred,udata,Fbp,F,fig);
% Sp=GetPrior1(fitP,stdPa,F.g2.L);




%.........................................................................
%model 2a
%Fit data's mean and std with std of llh(std) a hyperbolic function of
%coherence.
%4 free parameters for hyperbolic function,std of priors and motor noise
%.........................................................................
% display('Now fitting the model....')
% You are free to set (or not,[]) 9 free parameters.
% fitP=[];%[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
% [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut]=LSft2a(datafitm,fitP);
% display('Now drawing the predictions....')
% drawPred5(pred,udata,sdata,Fbp,F,fig)
% Sp=GetPrior2a(fitP,F.g2.L);%Get prior's representation



%.........................................................................
%model 3
%Fit data's mean with width of llh(std) as a free parameter and assuming
%the true priors
%.........................................................................
% [pred,udata,Fbp,fitP,R2]=LSft3(datafitm);




%.........................................................................
%model 4
%Fit data's mean and std with width of llh(std) as a free parameter and
%assuming the true priors
%.........................................................................
% [pred,udata,sdata,Fbp,fitP,R2]=LSft4(datafitm);
% [fig1]=drawPred4(pred,udata,sdata,Fbp,F,fig);




%.........................................................................
%model 5 (~30 min)
%Description:
%Fit data's mean and std with
%std of llh,std of priors and motor noise as 8 free parameters
%model 6 is fitted within the range [0 to 1000] because it can make
%prediction within this full range thanks to the analytical derivation.
%.........................................................................
fprintf(' \n %s \n','Now fitting the mean and std of the data ...')
%You are free to set (or not,[]) your own free parameters (row vector).
fitP=[];%[6 12 24 80 40 20 10 0.01]';
[logl,fitP,pred,data,d,coh,pstd,stdPa]=MLft5(fitP);
fprintf(' \n %s \n','Now drawing model predictions....')
[meansD,meansP,stdsD,stdsP]=drawPred5(pred,data,d,coh,pstd);
Sp=GetPrior5(fitP.p,stdPa,[0.74559 2.77 8.74 33.25],[.24 .12 .06]);
%.........................................................................
%model 6 (~17h not complete)
%Description: "I kinda see or I don't". Sometimes subjects kinda see
%(gaussian likelihood), sometimes they don't (flat likelihood). The
%proportion of each is implemented by a "fraction" parameter. It is fixed
%across directions coherences and priors.
%Fit mean and std with Bayesian inference
%std llh,std priors and as 7 free parameters
%there is no motor noise because I don't know yet how to implement it.
%If the code is ok data should follow the qualitative trends predicted with:
%fitP.p=[1/.24 1/.12 1/.06 80 40 20 10 5]';
%model 6 is fitted within the range [0.2 to 1000] because it can only make
%prediction within this range. Otherwise it outputs "NaN";
%Current problem is: how to calculate the variance of the posterior? The
%posterior is a bimodal distribution that is sometimes asymetric. I have
%noticed that when a distribution is asymetric the usual formula for
%variance gives a poor approximation of the true (population) variance.

%.........................................................................
% fprintf(' \n %s \n','Now fitting the model with the mean and std of the data ...')
% %You are free to set (or not,[]) the free parameters.
% fitP.p=[];%[1/.24 1/.12 1/.06 80 40 20 10 5]';
% [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa]=LSft6(datafitm,fitP);
% fprintf('\n %s \n','Now drawing model predictions....')
% fig1=drawPred6(pred,F,fig,datafitt);
% fprintf('\n %s \n','Now drawing prior....')
% Sp=GetPrior5(fitP.p,stdPa,F.g2.L,F.g1.L);


%full pdf
%[0.001 0.00100000000000011 0.001 241663.491210938 196214.551757812 183598.189453125 0.001 2873.797975]';
%R2=0.82

%llh's std fails-why?

%equation
%[16.0694875597013  22.739550821384 43.8654179979805 44.2595870883993 43.9075425556213 25.2721154896555 8.88794988907877 0.001]
%R2=0.922020837557773
%stdPa=[3.84426118499344 35.3881812184615 108.231155296127 85.9092739777982 59.5655105531189 44.3348025111216 26.2161732609227 30.1753390168421];



%Databank
function Makedatabank(data,FofLl2del,Ll2del)
global databank

%set the directory to analyse
pathFold=pwd;
d=dir(pathFold);

%collect subjects
datatmp.subjects=data.subjects;

%find the subfolders to analyse
for ii=1:numel(datatmp.subjects)
    isubType(ii)=find(strcmp({d.name},datatmp.subjects(ii)));
end

%get the names of the subfolders
Folds.name={d(isubType).name}';

%check if correct
disp(Folds.name)
%uiwait(msgbox('Check out the command window for the folder names.'));

%count the number of subfolders
Folds.nb=numel(Folds.name);

databank.data=[];

%loop over the subfolders to analyse
for j=1:Folds.nb
    
    %directory of the subfolder to open
    Folds.path(j)=strcat(pathFold,'/',Folds.name(j));
    
    %switch to the subfolder to open
    cd(Folds.path{j})
    
    %delete svn files (if exist)
    !rm ._*;
    
    %specify the files to load
    datadir=dir('*data*.mat');%directory
    
    %check if there are data in the directory
    if isempty(datadir)
        disp(strcat('No data were found in directory: ',Folds.name(j)))
        return
    end
    
    %collect the name of each file
    datatmp.nm={datadir.name};        %filenms
    
    %count the number of file
    datatmp.nb=numel(datatmp.nm);     %num`ber
    
    %set the variables in the databank (a column each)
    filedetails={};
    Trials={};
    session={};
    runr={};
    es_coor={};
    Pstd={};
    priormean={};
    sample_dir={};
    coh={};
    
    
    %check if data have been specified in the function argin
    if ~isempty(datatmp.nm)
        datalisting=datatmp.nm;
        %    datalisting={datalisting};
        %    display(["--- A databank is being created with the data specified ----']);
        
        %if data have not been specified,gather data from directory
    else
        %remove possible svn files
        datalisting=dir('*data*.mat');
        datalisting={datalisting.name};
        %    display(['--- No data have been specified. A databank is being created with all data in the directory and analysed ---']);
    end
    
    %tic
    
    %loop over the files and collect their data
    for i=1:numel(datalisting)
        
        %load the files and get data from 'task.mat' in the workspace
        load(datalisting{i});
        datai=getTaskParameters(myscreen,task);%speed consuming
        
        %get the estimated cartesian coordinates(data)
        es_coor_i=datai{2}.randVars.prodcoor';
        es_coor=[es_coor;es_coor_i];
        
        %convert coordinates to angles (degree)
        coortmp=cell2mat(es_coor);
        esf4Visu=num2cell(getangle(coortmp(:,1),coortmp(:,2)));
        
        %calculate the number of trials
        numTrials=numel(es_coor_i);
        Trials_i=num2cell(1:numTrials)';
        Trials=[Trials;Trials_i];
        
        %get the file name
        filedetails_i=repmat(datalisting(i),numTrials,1);
        filedetails=[filedetails;filedetails_i];
        
        %get the run
        if isempty(strfind(filedetails_i{1},'run'))
            disp('the filename does not contains "run" information')
            return
        end
        run_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'run')+3:strfind(filedetails_i{1},'run')+4));%run
        runr=[runr;num2cell(repmat(run_thisT,numTrials,1))];
        
        %get the session
        if isempty(strfind(filedetails_i{1},'sess'))
            disp(['--- the filename does not contains "session" information ---'])
            return
        end
        session_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'sess')+4));%session
        session=[session;num2cell(repmat(session_thisT,numTrials,1))];
        
        %get the std of the prior from the file nm;
        if isempty(strfind(filedetails_i{1},'Pstd'))
            disp(['--- the filename does not contains "Pstd" information ---'])
            return
        end
        Pstd_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2));%std of the prior
        
        %get the std of the prior from the file nm when prior is uniform;
        if isempty(strfind(filedetails_i{1},'inf'))==0
            Pstd_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'Pstd')+numel('Pstd'):strfind(filedetails_i{1},'Pstd')+numel('Pstd')+2));%std of the prior
        end
        
        %check if the Prior's std was found in the name
        if isnan(Pstd_thisT)
            disp(['--- the code does not find the value of "Pstd" in the filename ---'])
            return
        end
        Pstd=[Pstd;num2cell(repmat(Pstd_thisT,numTrials,1))];
        
        %get the prior's mean
        if isempty(strfind(filedetails_i{1},'mean'))
            disp(['--- the filename does not contains "Prior mean" information ---'])
            return
        end
        priormean_thisT=str2double(filedetails_i{1}(strfind(filedetails_i{1},'mean')+4:strfind(filedetails_i{1},'mean')+6));%run
        priormean=[priormean;num2cell(repmat(priormean_thisT,numTrials,1))];
        
        %get the Displayed directions (group 1)
        %when design is balanced across coherences (same trial number).
        if isfield(datai{1,2}.parameter,'dir')
            sample_dir_i=num2cell(datai{1,2}.parameter.dir)';
            sample_dir=[sample_dir;sample_dir_i];
            %when design is adjusted according to coherences (different trial nb).
        elseif isfield(datai{1,2}.randVars,'myRandomDir')
            sample_dir_i=num2cell(datai{1,2}.randVars.myRandomDir)';
            sample_dir=[sample_dir;sample_dir_i];
        end
        %get the sample coherences (group 2)
        %when design is balanced across coherences (same trial number).
        if isfield(datai{1,2}.parameter,'coherence')
            coh_i=num2cell(datai{1,2}.parameter.coherence)';
            coh=[coh;coh_i];
            %when design is adjusted according to coherences (different trial nb).
        elseif isfield(datai{1,2}.randVars,'myRandomCoh')
            coh_i=num2cell(datai{1,2}.randVars.myRandomCoh)';
            coh=[coh;coh_i];
        end
    end
    %toc
    %create our databank
    %enter the names of the variable in the columns
    databank.nm=[
        {'filedetails'  },...
        {'run'          },...
        {'session'      },...
        {'Trials'       },...
        {'es_coor'},...
        {'Pstd'         },...
        ('priormean'    ),...
        {'sample_dir'   },...
        {'coh'          },...
        {'est_dir' }];
    
    %store the data in each column
    databanki.data=[
        filedetails(:)  ...
        runr(:)         ...
        session(:)      ...
        Trials(:)       ...
        es_coor(:)     ...
        Pstd(:)         ...
        priormean(:)    ...
        sample_dir(:)   ...
        coh(:)          ...
        esf4Visu];
    
    %store data over subjects
    databank.data=[databank.data;databanki.data];
end

%data are organized and saved in a file called datbank in the directory
%switch back to the mother directory
cd(pathFold)

%discard a particular condition from the databank
FofLl2deli=strcmp(databank.nm,{FofLl2del});
Ll2deli=cell2mat(databank.data(:,FofLl2deli))==Ll2del;
databank.data(Ll2deli,FofLl2deli)={[]};

%backup the databank
save('datbank','databank');


%Sort factors
function [datafitm,datafitt,F]=initFs(Fs,fig)
global databank

%data
es_data=cell2mat(databank.data(:,(strcmp(databank.nm,'es_coor'))==1));
es.dir=cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1));

%factors
%F1
F.g1.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,Fs{1}  ))==1));
F.g1.nm=Fs{1};
F.g1.L=unique(F.g1.thisT);
F.g1.L=sort(F.g1.L,'descend');
F.g1.numL=numel(F.g1.L);
clear i
for i=1:F.g1.numL
    F.g1.posLli(i)={find(F.g1.thisT==F.g1.L(i))};
end

%F2
F.g2.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,Fs{2}  ))==1));
F.g2.nm=Fs{2};
F.g2.L=unique(F.g2.thisT);
F.g2.L=sort(F.g2.L,'descend');
F.g2.numL=numel(F.g2.L);
for i=1:F.g2.numL
    F.g2.posLli(i)={find(F.g2.thisT==F.g2.L(i))};
end

%F3
F.g3.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,Fs{3}  ))==1));
F.g3.nm=Fs{3};F.g3.L=unique(F.g3.thisT);
F.g3.L=sort(F.g3.L,'ascend');
F.g3.numL=numel(F.g3.L);
for i=1:F.g3.numL
    F.g3.posLli(i)={find(F.g3.thisT==F.g3.L(i))};
end

%position of conditions
for k=1:F.g1.numL
    for j=1:F.g2.numL
        for i=1:F.g3.numL
            F.g1g2g3.posLli(i,j,k)=...
                {intersect( ...
                intersect(F.g1.posLli{k},F.g2.posLli{j}),...
                F.g3.posLli{i})};
        end
    end
end

[datafitm,datafitt,F]=processData1(es,es_data,F);
%[datafitm,datafitt,F]=processData2(es_data,F);


%Process data
function [datafitm,datafitt,F]=processData1(es,es_data,F)
%Calculate coordinates of average estimated directions for each condition
%organize data in following order: group 1(subplots) - group 2(colors) -
%group 3(x-axis). Each cell contains the repetitions of a condition.
for k=1:F.g1.numL
    for j=1:F.g2.numL
        for i=1:F.g3.numL
            
            %calculate mean & std
            es_sta{i,j,k}=statcircular(es_data(F.g1g2g3.posLli{i,j,k},:));
            
            %get mean estimates (degree)
            es.mean(i,j,k)=es_sta{i,j,k}.deg.mean;
            es.std(i,j,k)=es_sta{i,j,k}.deg.std;
            
            %get mean coordinates
            es.meanCo{i,j,k}=es_sta{i,j,k}.coord.mean;
            
            %get factors & levels
            es.Fbp1(i,j,k)=F.g1.L(k);
            es.Fbp2(i,j,k)=F.g2.L(j);
            es.Fbp3(i,j,k)=F.g3.L(i);
            
            %get motion direction coordinates
            %radius
            r=2.5;
            disp.coord{i,j,k}=polar2cartesian(es.Fbp3(i,j,k),r);
        end
    end
end

%mean data sorted by factors
datafitm.data=[es.mean(:),...
    es.std(:),...
    es.Fbp1(:),...
    es.Fbp2(:),...
    es.Fbp3(:)];

datafitm.nm={'mean','std','F1','F2','F3'};

%trial-data
datafitt.data=[{es.dir},...
    {es.mean},...
    {es.std},...
    {es.Fbp1},...
    {es.Fbp2},...
    {es.Fbp3}];
datafitt.nm={'dir','mean','std','F1','F2','F3'};
function [datafitm,datafitt,F]=processData2(es_data,F)
%This function calculates the estimates' average distance of the prior.
%It is not clear what this variable means or what assumptions must be made
%to work with it so I'll avoid it.
%Calculating estimates' average distances to the prior is different from
%measuring the distance of estimates' vector average to the prior and so,
%yield different results. What this difference means is not clear.

%radius
r=2.5;

%prior's coordinates
prior.coo=polar2cartesian(225,r);

%motion direction coordinates
%sample-wise
F.g3.Lcoo=polar2cartesian(F.g3.L,r);
%trial-wise
F.g3.thisTcoo=polar2cartesian(F.g3.thisT,r);

%express data in distance to the prior
for i=1:size(es_data,1)
    %estimated
    es.dist(i)=vectors2signedAngle(es_data(i,:),prior.coo);
    %motion direction
    F.g3.distthisT(i)=vectors2signedAngle(F.g3.thisTcoo(i,:),prior.coo);
end
%motion direction
for i=1:size(F.g3.Lcoo,1)
    F.g3.distL(i)=vectors2signedAngle(F.g3.Lcoo(i,:),prior.coo);
end

%mean,std,g1,g2,g3
for k=1:F.g1.numL
    for j=1:F.g2.numL
        for i=1:F.g3.numL
            es.mean(i,j,k)=mean(es.dist(F.g1g2g3.posLli{i,j,k}));
            es.std(i,j,k)=std(es.dist(F.g1g2g3.posLli{i,j,k}));
            es.Fbp1(i,j,k)=F.g1.L(k);
            es.Fbp2(i,j,k)=F.g2.L(j);
            es.Fbp3(i,j,k)=F.g3.distL(i);
        end
    end
end

%data for fit to mean
datafitm.data=[es.mean(:),...
    es.std(:),...
    es.Fbp1(:),...
    es.Fbp2(:),...
    es.Fbp3(:)];
datafitm.nm={'mean','std','F1','F2','F3'};

%data for trial-to-trial fit
datafitt=[{es.dist'},...
    {es.mean},...
    {es.std},...
    {es.Fbp1},...
    {es.Fbp2},...
    {es.Fbp3}];
datafitt.nm={'distance','mean','std','F1','F2','F3'};


%Fit
function [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa]=LSft(datafitm,fitP)
%Remove "NaN" data
l=datafitm.data(:,1);
datafitm.data(arrayfun(@(l) any(isnan(l)),l),:)=[];

%data
udata=datafitm.data(:,1);
sdata=datafitm.data(:,2);

%Factors
%g1
F.g1=datafitm.data(:,3);
F.g1L=unique(F.g1);
F.g1L=sort(F.g1L,'descend');

%g2
F.g2=datafitm.data(:,4);
F.g2L=unique(F.g2);
F.g2L=sort(F.g2L,'descend');

%g3
F.g3=datafitm.data(:,5);
F.g3L=unique(F.g3);
F.g3L=sort(F.g3L,'ascend');

%Store
Fbp=[F.g1 F.g2 F.g3];

%options
options=optimset('Display','off','MaxFunEvals',50000,...
    'TolFun',1e-10,...
    'MaxIter',50000,...
    'Algorithm','active-set');

%If there is no parameters,fit the model to the data
if isempty(fitP)==1
    
    %Initial parameters
    k0_p1=0:40:119;%before 0.05:0.05:0.15;
    k0_p2=0:40:119;
    k0_p3=0:40:119;
    k0_p4=0:40:119;
    
    %Fit
    %iteration=0;
    tic
    for i=1:numel(k0_p1)
        for j=1:numel(k0_p2)
            for k=1:numel(k0_p3)
                for l=1:numel(k0_p4)
                    
                    %Fit
                    [fitPtmp,~,~,~,~,~,Hessian]=fmincon( @(fitPtmp) makeSSE1(udata,Fbp,...
                        fitPtmp),...
                        [k0_p1(i);k0_p2(j);k0_p3(k);k0_p4(l)],...
                        [],[],[],[],...
                        [0 0 0 0],...
                        [inf inf inf inf],...
                        [],...
                        options);
                    
                    %Store
                    fitPbkp{i,j,k,l}=fitPtmp;
                    
                    %Calculate SSE
                    SSE_bkp(i,j,k,l)=makeSSE1(udata,Fbp,fitPtmp);
                    
                    toc
                    %    %Check
                    %    fprintf('%6.2f \n',SSE_bkp(i,j,k,l,m,n,e,f))
                    %    iteration=iteration + 1;
                    %    hold all
                    %    plot(iteration,SSE_bkp(i,j,k,l,m,n,e,f),'o','markerfacecolor','b',...
                    %        'markersize',13)
                    %    drawnow
                end
            end
        end
    end
    
    %Get the lower SSE
    [minSSE,position]=min(SSE_bkp(:));
    [i,j,k,l]=ind2sub(size(SSE_bkp),position);
    
    %Get the best parameters
    fitP=fitPbkp{i,j,k,l};
    
    %Get the model's predictions
    [pred,Factors]=makePre1(F,fitP);
    
    %data-based model predictions
    [~,pred00]=makeSSE1(udata,Fbp,fitP);
    
    %Get the R^2
    %%method 1
    R2=makeR2([udata;sdata],minSSE);
    %%method 2
    %[~,pred2]=makeSSE5(udata,sdata,Fbp,fitP);
    %R2=makeR22([udata;sdata],[pred2.mean';pred2.std']);
    
    %Store output
    dataOut=udata;
    predOut=pred00.mean;
    
    %std of model parameters
    stdPa=GetStdP(Hessian,dataOut,predOut,fitP);
    
    %If free parameters are input,only draw predictions
elseif isempty(fitP)==0
    
    %Calculate the SSE
    [~,pred00]=makeSSE1(udata,Fbp,fitP);
    
    %Calculate the R2
    R2=makeR22(udata,pred00.mean);
    
    %Store output
    dataOut=udata;
    predOut=pred00.mean';
    
    %Get model predictions on a larger space
    pred=makePre1(F,fitP);
    
    stdPa=[];
end
function [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut]=LSft2a(datafitm,fitP)
%GOAL
%Fit a model where the width of the llh is a hyperbolic function of coherence.
%Busse,L.,Ayaz,A.,Dhruv,N. T.,Katzner,S.,Saleem,A. B.,SchoLinck,M. L.,et al. (2011).
%The Detection of Visual Contrast in the Behaving Mouse. Journal of Neuroscience,31(31),11351?11361.
%doi:10.1523/JNEUROSCI.6689-10.2011
%Albrecht,D. G.,& Hamilton,D. B. (1982). Striate cortex of monkey and cat: contrast response
%function. J Neurophysiol.
%Least Square Optimization (LSO)

%INPUTS
%"datafitm": matrix of data and factors
%col1: data
%col2: factor 1
%col3: factor 2
%col4: factor 3
%fitP: vector of free parameters or empty matrix
%e.g.,fitP=[0.1 0.2 0.3 0.4];
%e.g.,fitP=[];


%Get data & factors
%.........................................................................................................................--

%Remove "NaN" data
l=datafitm.data(:,1);
datafitm.data(arrayfun(@(l) any(isnan(l)),l),:)=[];

%Collect data
udata=datafitm.data(:,1);
sdata=datafitm.data(:,2);

%Factors
%g1
F.g1=datafitm.data(:,3);
F.g1L=unique(F.g1);
F.g1L=sort(F.g1L,'descend');%order

%g2
F.g2=datafitm.data(:,4);
F.g2L=unique(F.g2);
F.g2L=sort(F.g2L,'descend');%order

%g3
F.g3=datafitm.data(:,5);
F.g3L=unique(F.g3);
F.g3L=sort(F.g3L,'ascend');%order

%Store
Fbp=[F.g1 F.g2 F.g3];


%options
options=optimset('Display','off','MaxFunEvals',50000,...
    'TolFun',1e-06,...
    'MaxIter',50000,...
    'Algorithm','active-set');

%If there is no parameters,fit the model to the data
if isempty(fitP)==1
    
    %Initial parameters
    R0_0=0%:0.05:0.1;
    Rmax_0=1%:0.05:0.1;
    n_0=1%:0.05:0.1;
    C50_0=1%:0.05:0.1;
    sp0_1=1:500:1000;
    sp0_2=1:500:1000;
    sp0_3=1:500:1000;
    sp0_4=1:500:1000;
    sM0=0%:0.05:0.1;
    
    %Fit
    %.........................................................................................................................--
    %iteration=0;
    %Loop over free parameters
    for i=1:numel(R0_0)%e.g.,
        for j=1:numel(Rmax_0)%e.g.,
            for k=1:numel(n_0)%e.g.,motor noise
                for l=1:numel(C50_0)%e.g.,
                    for m=1:numel(sp0_1)%e.g.,prior 1
                        for n=1:numel(sp0_2)%e.g.,prior 2
                            for o=1:numel(sp0_3)%e.g.,prior 3
                                for p=1:numel(sp0_4)%e.g.,prior 4
                                    for q=1:numel(sM0)%e.g.,motor noise
                                        
                                        %Fit
                                        fitPtmp=fmincon( @(fitPtmp) makeSSE2a(udata,sdata,Fbp,...
                                            fitPtmp),...
                                            [R0_0(i);Rmax_0(j);n_0(k);C50_0(l);sp0_1(m);sp0_2(n);sp0_3(o);sp0_4(p);sM0(q)],...
                                            [],[],[],[],...
                                            [.1 .1 .1 .1 .1 .1 .1 .1 .1],...
                                            [inf inf inf inf inf inf inf inf inf],...
                                            [],...
                                            options);
                                        
                                        %Store
                                        fitPbkp{i,j,k,l,m,n,o,p,q}=fitPtmp;
                                        
                                        %Calculate SSE
                                        SSE_bkp(i,j,k,l,m,n,o,p,q)=makeSSE2a(udata,sdata,Fbp,fitPtmp);
                                        
                                        %    %Check
                                        %    fprintf('%6.2f \n',SSE_bkp(i,j,k,l,m,n,e,f))
                                        %    iteration=iteration + 1;
                                        %    hold all
                                        %    plot(iteration,SSE_bkp(i,j,k,l,m,n,e,f),'o','markerfacecolor','b',...
                                        %        'markersize',13)
                                        %    drawnow
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %Get the lower SSE
    [minSSE,position]=min(SSE_bkp(:));
    [i,j,k,l,m,n,o,p,q]=ind2sub(size(SSE_bkp),position);
    
    %Get the best parameters
    fitP=fitPbkp{i,j,k,l,m,n,o,p,q};
    
    %Get the model's predictions
    [pred,~]=makePre2a(F,fitP);
    
    %Get the R^2
    %%method 1
    R2=makeR2([udata;sdata],minSSE);
    %%method 2
    %[~,pred2]=makeSSE5(udata,sdata,Fbp,fitP);
    %R2=makeR22([udata;sdata],[pred2.mean';pred2.std']);
    
    %Store output
    dataOut=udata;
    predOut=pred.mean;
    
    %If free parameters are input,only draw predictions
elseif isempty(fitP)==0
    
    %Calculate the SSE
    [~,pred00]=makeSSE1(udata,Fbp,fitP);
    
    %Calculate the R2
    R2=makeR22(udata,pred00.mean);
    
    %Store output
    dataOut=udata;
    predOut=pred00.mean';
    
    %Get model predictions on a larger space
    pred=makePre2a(F,fitP);
    
end
function [pred,udata,Fbp,fitP,R2,fig1]=LSft3(datafitm)

%Fit data for each prior separately.
%Parameters search (model fitting by Least Square Optimization (LSO))

%e.g.,of the parameters used:
%c=0.12;
%s=2;
%sPr=20;
%uPr=225;
%uLl=155:10:295;
%%two free parameters:
% sl,sp


%Remove "NaN" data
l=datafitm.data(:,1);
datafitm.data(arrayfun(@(l) any(isnan(l)),...
    l),:)=[];


%Fit sl for each coherence
%for data(coh==1),set the priors as true and fit sl.

%Collect the data for coherence ii
Ftofit.labels=datafitm.data(:,3);%(e.g.,coh)
Ftofit.L=sort(unique(Ftofit.labels),'descend');

%Initial parameters
k0=0.05:0.05:20;


fitPBkp=[];
%Loop over coherences
for ii=1:numel(Ftofit.L);
    
    %.........................................................................................................................--
    %Select variables
    %.........................................................................................................................--
    %Data
    udata{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii),1);
    
    %Conditions
    %g1
    F.g1{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii) ,3);
    %g2
    F.g2{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii) ,4);
    %g3
    F.g3{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii) ,5);
    %Store
    Fbp{ii}=[F.g1{ii}  F.g2{ii}  F.g3{ii}];
    
    %.........................................................................................................................--
    %Fit data
    %.........................................................................................................................--
    
    %options
    options=optimset('Display','off','MaxFunEvals',50000,...
        'TolFun',1e-06,...
        'MaxIter',50000,...
        'Algorithm','active-set');
    
    %Loop over the initial parameters
    for i=1:numel(k0)
        
        %Fitting
        fitPtmp=fmincon( @(fitPtmp) makeSSE3(udata{ii},Fbp{ii},...
            fitPtmp),...
            k0(i),...
            [],[],[],[],0,1000,[],...
            options);
        
        %Store free parameters
        fitPBkup{1,ii}(i)=fitPtmp;
        
        %Calculate the SSE
        SSE_bkp{ii}(i)=makeSSE3(udata{ii},Fbp{ii},fitPBkup{1,ii}(i));
        
    end
    
    %.........................................................................................................................--
    %Draw SSE,predictions and data
    %.........................................................................................................................--
    %Get the lower SSE
    minSSE{ii}=min(SSE_bkp{ii});
    [~,idxcol]=min(SSE_bkp{ii});
    
    %Get the parameter for the lower SSE and the factor
    fitP{ii,1}=fitPBkup{ii}(idxcol);
    fitP{ii,2}=Ftofit.L(ii);
    
    %Compute the model's best predictions
    [~,pred{ii},~]=makeSSE3(udata{ii},Fbp{ii},fitP{ii,1});
    
    %Calculate R^2,the percent of variance in the data explained by the
    %model.
    R2{ii}=makeR2(udata{ii},minSSE{ii});
    
    
    %Check the convergence of the optimization
    figure('color',[1 1 1]);
    subplot(211);
    title('SSE','fontsize',12)
    hold all
    plot(k0,SSE_bkp{ii},'k');
    plot(k0(idxcol),minSSE{ii},'O','markerfacecolor','r','markersize',10)
    subplot(212);
    title('sl','fontsize',12)
    hold all
    plot(k0,fitPBkup{1,ii},'k');
    plot(k0(idxcol),fitP{ii,1},'O','markerfacecolor','r','markersize',10)
    
end
function [pred,udata,sdata,Fbp,fitP,R2,fig1]=LSft4(datafitm)
global fig
%Fit mean and variance with std(llh) as a free parameter and assuming the true priors
%Parameters search (model fitting by Least Square Optimization (LSO))

%e.g.,of the parameters used:
%c=0.12;
%s=2;
%sPr=20;
%uPr=225;
%uLl=155:10:295;
%%two free parameters:
% sl,sp


%Remove "NaN" data
l=datafitm.data(:,1);
datafitm.data(arrayfun(@(l) any(isnan(l)),...
    l),:)=[];

%Fit sl for each coherence
%for data(coh==1),set the priors as true and fit sl.

%Collect the data for coherence ii
Ftofit.labels=datafitm.data(:,3);%(e.g.,coh)
Ftofit.L=sort(unique(Ftofit.labels),'descend');

%Initial parameters
k0=0.05;%:0.05:20;


fitPBkp=[];
%Loop over coherences
for ii=1:numel(Ftofit.L);
    
    %.........................................................................................................................--
    %Select variables
    %.........................................................................................................................--
    %Data
    udata{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii),1);
    sdata{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii),2);
    
    %Conditions
    %g1
    F.g1{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii) ,3);%e.g.,coh
    %g2
    F.g2{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii) ,4);
    %g3
    F.g3{ii}=datafitm.data(Ftofit.labels==Ftofit.L(ii) ,5);
    %Store
    Fbp{ii}=[F.g1{ii}  F.g2{ii}  F.g3{ii}];
    
    %.........................................................................................................................--
    %Fit data
    %.........................................................................................................................--
    
    %options
    options=optimset('Display','off','MaxFunEvals',50000,...
        'TolFun',1e-06,...
        'MaxIter',50000,...
        'Algorithm','active-set');
    %Loop over the initial parameters
    for i=1:numel(k0)
        
        %Fitting
        fitPtmp=fmincon( @(fitPtmp) makeSSE4(udata{ii},sdata{ii},Fbp{ii},...
            fitPtmp),...
            k0(i),...
            [],[],[],[],0,1000,[],...
            options);
        
        %Store free parameters
        fitPBkup{1,ii}(i)=fitPtmp;
        
        %Calculate the SSE
        SSE_bkp{ii}(i)=makeSSE4(udata{ii},sdata{ii},Fbp{ii},fitPBkup{1,ii}(i));
        
    end
    
    %.........................................................................................................................--
    %Draw SSE,predictions and data
    %.........................................................................................................................--
    %Get the lower SSE
    minSSE{ii}=min(SSE_bkp{ii});
    [~,idxcol]=min(SSE_bkp{ii});
    
    %Get the parameter for the lower SSE and the factor
    fitP{ii,1}=fitPBkup{ii}(idxcol);
    fitP{ii,2}=Ftofit.L(ii);
    
    %Compute the model's best predictions
    [~,pred{ii},~]=makeSSE4(udata{ii},sdata{ii},Fbp{ii},fitP{ii,1});
    
    %Calculate R^2,the percent of variance in the data explained by the
    %model.
    R2{ii}=makeR2([udata{ii};sdata{ii}],minSSE{ii});
    
    %Check the convergence of the optimization
    figure('color',[1 1 1]);
    subplot(211);
    title('SSE','fontsize',12)
    hold all
    plot(k0,SSE_bkp{ii},'k');
    plot(k0(idxcol),minSSE{ii},'O','markerfacecolor','r','markersize',10)
    subplot(212);
    title('sl','fontsize',12)
    hold all
    plot(k0,fitPBkup{1,ii},'k');
    plot(k0(idxcol),fitP{ii,1},'O','markerfacecolor','r','markersize',10)
    xlabel('initial values')
end
function [logl,fitP,pred,data,disp,coh,pstd,stdPa]=MLft5(fitP)
global databank

%data
% estimated
data=cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1));
data=round(data);
i=isnan(data)==0;
data=data(i);

%Fixed parameters
%direction
disp=cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));
disp=disp(i);

%coherence
coh=cell2mat(databank.data(:,(strcmp(databank.nm,'coh'))==1));
coh=coh(i);

%prior std
pstd=cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
pstd=pstd(i);

% % %................
%test (to remove)
% load data
% data=upo;
% disp=d';
% coh=coh';
% pstd=pstd';
% %.................

%options
%search with fmincon algorithm interior-point
options=optimset('Display','iter',...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'Algorithm','interior-point');

%Search with fminsearch Nelder-Mead
% Options=optimset('Display','iter','MaxFunEvals',1000,...
%     'TolFun',1e-10,...
%     'TolX',1e-10,...
%     'MaxIter',5000);

%Search with fminsearch LevenbergMarquardt (doesn't work well)
% Options=optimset('Display','iter','MaxFunEvals',5000,...
%     'TolFun',1e-04,...
%     'TolX',1e-06,...
%     'MaxIter',5000,...
%     'Algorithm','levenberg-marquardt',...
%     'ScaleProblem','Jacobian');

tic
%If there are no parameters, fit.
if isempty(fitP)==1
    
    %Initial parameters: found by manually matching raw data with
    %simulations predictions.
    sl1_0=19;%[250 350 450];%350;%k(1);
    sl2_0=9;%[14 24 34];%24;%k(2);
    sl3_0=1.7;%[1.5 7.5 17.5];%7.5;%k(3);
    sp1_0=1.05;%[1.3 4.3 14.3];%4.3;%k(4);
    sp2_0=1.7;%[1.9 8.9 18.9];%8.9;%k(5);
    sp3_0=2.5;%[4.5 14.5 24.5];%14.5;%k(6);
    sp4_0=7;%[23 33 43];%33;%k(7);
    Pm=0.005;
    
    
    %Fitting
    %k must be in the range [0:709]. when k>709, von mises probabilities
    %are inf.
    %k must be in the range [0:700]. when k>700, besseli is inf amd mPdfs
    %is 0. When normalizing the von mises to probabilities we get NaN
    %values(0/sum(0)).
    %k must be >0 because when both k are 0, Bayesian equation produces NaN
    %in atan2 (kl/kp).
    %So k is between 0.001 and 700.
    tic;
    fitPbkp=nan(8,numel(sl1_0),numel(sl2_0),numel(sl3_0),numel(sp1_0),numel(sp2_0),numel(sp3_0),numel(sp4_0),numel(Pm));
    logL_bkp=nan(numel(sl1_0),numel(sl2_0),numel(sl3_0),numel(sp1_0),numel(sp2_0),numel(sp3_0),numel(sp4_0),numel(Pm));
    for i=1:numel(sl1_0)
        for j=1:numel(sl2_0)
            for k=1:numel(sl3_0)
                for l=1:numel(sp1_0)
                    for m=1:numel(sp2_0)
                        for n=1:numel(sp3_0)
                            for e=1:numel(sp4_0)
                                for f=1:numel(Pm)

                                %Fit (interior-point)
                                [fitPtmp,logL,exitflag,output,~,~,Hessian]=fmincon( @(fitPtmp) makeSSE5(data,disp,coh,pstd,...
                                    fitPtmp),...
                                    [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e);Pm(f)],...
                                    [],[],[],[],...
                                    [0 0 0 0 0 0 0 0],...
                                    [50 50 50 50 50 50 50 1],...
                                    [],...
                                    options);
                                
                                %search way 2 (Nelder-Mead)
                                %[fitPtmp,logL,exitflag,output]=fminsearch(@(fitPtmp) makeSSE5(data,disp,coh,pstd,...
                                %    fitPtmp),...
                                %    [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e)],...
                                %    Options)
                                
                                %search way 3 (Levenberg-Marquardt)
                                %[fitPtmp,~,logL,exitflag,output,~,jacobian]=lsqnonlin(@(fitPtmp) makeSSE5(data,disp,coh,pstd,...
                                %     fitPtmp),...
                                %     [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e)],...
                                %     [],...
                                %     [],...
                                %    Options)
                                
                                %Fit parameters and SSE
                                fitPbkp(:,i,j,k,l,m,n,e,f)=fitPtmp;
                                logL_bkp(i,j,k,l,m,n,e,f)=logL;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %Max likelihood. To have an idea of what log likelihood we should
    %obtain our reasoning is that if one cannot predict at all motion
    %direction at a given trial, the probability that any given direction
    %be chosen is 1/360. (we consider estimate with 1 degree resolution).
    %So loglikelihood should be log(1/360)=-5.88. This mutiplied (because
    %of log) by our ~5500 trials, we get logL~32000. Any ability to
    %somewhat predict better motion direction better than chance yields
    %soemthing lower than 32000.
    [logl,position]=min(logL_bkp(:));
    [i,j,k,l,m,n,e,f]=ind2sub(size(logL_bkp),position);
    
    %best parameters
    fitP.p=fitPbkp(:,i,j,k,l,m,n,e,f);
    fitP.nm={'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','Pmotor'};
    
    %general predictions
    pred=makePre5(disp,coh,pstd,fitP.p);
    
    %std of model parameters
    stdPa=GetStdP(Hessian,data,pred,fitP.p);
end
function [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa]=LSft6(datafitm,fitP)

%Remove missing
l=datafitm.data(:,1);
datafitm.data(arrayfun(@(l) any(isnan(l)),l),:)=[];

%data
udata=datafitm.data(:,1);
sdata=datafitm.data(:,2);

%Conditions
%g1
F.g1=datafitm.data(:,3);
F.g1L=unique(F.g1);
F.g1L=sort(F.g1L,'descend');

%g2
F.g2=datafitm.data(:,4);
F.g2L=unique(F.g2);
F.g2L=sort(F.g2L,'descend');

%g3
F.g3=datafitm.data(:,5);
F.g3L=unique(F.g3);
F.g3L=sort(F.g3L,'ascend');

%Store
Fbp=[F.g1 F.g2 F.g3];

%options
options=optimset('Display','off','MaxFunEvals',10000,...
    'TolFun',1e-06,...
    'MaxIter',10000,...
    'Algorithm','active-set');

%If there is no parameters,fit the model to the data
if isempty(fitP.p)==1
    
    %Initial parameters
    sl1_0=1:100;%:164;
    sl2_0=1:100;%:164;
    sl3_0=1:100;%:164;
    sp1_0=1:100;%:164;
    sp2_0=1:100;%:164;
    sp3_0=1:100;%:164;
    sp4_0=1:100;%:164;
    sM_0=1;
    
    %.......
    %Fitting
    %.......
    %iteration=0;
    tic
    for i=1:numel(sl1_0)
        for j=1:numel(sl2_0)
            for k=1:numel(sl3_0)
                for l=1:numel(sp1_0)
                    for m=1:numel(sp2_0)
                        for n=1:numel(sp3_0)
                            for e=1:numel(sp4_0)
                                for f=1:numel(sM_0)
                                    
                                    %Fit
                                    %0.2 to 1000 is chosen because model 6
                                    %can only make prediction within this
                                    %range. Otherwise it outputs "NaN";
                                    [fitPtmp,~,~,~,~,~,Hessian]=fmincon( @(fitPtmp) makeSSE6(udata,sdata,Fbp,...
                                        fitPtmp),...
                                        [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e);sM_0(f)],...
                                        [],[],[],[],...
                                        [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2],...
                                        [1000 1000 1000 1000 1000 1000 1000 1000],...
                                        [],...
                                        options);
                                    
                                    %fit parameters and SSE
                                    fitPbkp{i,j,k,l,m,n,e,f}=fitPtmp;
                                    SSE_bkp(i,j,k,l,m,n,e,f)=makeSSE6(udata,sdata,Fbp,fitPtmp);
                                    
                                    toc
                                    %                                    %Check
                                    %                                    fprintf('%6.2f \n',SSE_bkp(i,j,k,l,m,n,e,f))
                                    %                                    iteration=iteration + 1;
                                    %                                    hold all
                                    %                                    plot(iteration,SSE_bkp(i,j,k,l,m,n,e,f),'o','markerfacecolor','b',...
                                    %                                        'markersize',13)
                                    %                                    drawnow
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    %lower SSE
    [minSSE,position]=min(SSE_bkp(:));
    [i,j,k,l,m,n,e,f]=ind2sub(size(SSE_bkp),position);
    
    %best parameters
    fitP.p=fitPbkp{i,j,k,l,m,n,e,f};
    fitP.nm={'24','12','6','80','40','20','10','motor'};
    
    %general predictions
    [pred,~]=makePre(F,fitP.p);
    
    %data-based predictions
    [~,pred00]=makeSSE6(udata,sdata,Fbp,fitP.p);
    
    %R^2
    %%way 1
    R2=makeR2([udata;sdata],minSSE);
    %%way 2
    %[~,pred2]=makeSSE5(udata,sdata,Fbp,fitP);
    %R2=makeR22([udata;sdata],[pred2.mean';pred2.std']);
    
    %output
    dataOut=[udata;sdata];
    predOut=[pred00.mean;pred00.std];
    
    %std of model parameters
    stdPa=GetStdP(Hessian,dataOut,predOut,fitP.p);
    
    %If free parameters are input, only draw predictions
elseif isempty(fitP.p)==0
    
    %SSE
    [~,pred00]=makeSSE6(udata,sdata,Fbp,fitP.p);
    
    %R2
    R2=makeR22([udata;sdata],[pred00.mean;pred00.std]);
    
    %output
    dataOut=[udata;sdata];
    predOut=[pred00.mean';pred00.std'];
    
    %general predictions
    [pred,~]=makePre(F,fitP.p);
    
    %std of model parameters is not meaningful
    stdPa=[];
end


%SSE/logL
function [SSE,pred,fitP]=makeSSE1(udata,Fbp,fitP)

%Representations
%llh
%mean
uLl=Fbp(:,3);
%coherence
c=Fbp(:,1);

%prior
%mean
uPr=225;
%std
k=nan(numel(udata),1);
k(Fbp(:,2)==80)=fitP(1);
k(Fbp(:,2)==40)=fitP(2);
k(Fbp(:,2)==20)=fitP(3);
k(Fbp(:,2)==10)=fitP(4);

%Bayes infer
%https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
%mean
uPo=(1./(1+(k./c).^2)).*uLl+(1./(1+(c./k).^2)).*uPr;

%calculate the SSE and predictions
SSE=nansum((udata-uPo).^2);
pred.mean=uPo;
function [SSE,pred,fitP]=makeSSE2a(udata,sdata,Fbp,fitP)



%Representations
%...............
%Task parameters
c =Fbp(:,1);%coherence
spe=Fbp(:,2);%std of the prior
uPr=225;%Mean of the prior
uLl=Fbp(:,3);%Mean of the likelihood

%Free parameters
%llh
R0=fitP(1);
Rmax=fitP(2);
n=fitP(3);
C50=fitP(4);
%prior
sp1=fitP(5);%std of prior 1
sp2=fitP(6);
sp3=fitP(7);
sp4=fitP(8);
%motor
sM=fitP(9);%motor noise

%.............................................................................
%Equation when the noise of LLH is an hyperbolic function of the coherence
%Busse,L.,Ayaz,A.,Dhruv,N. T.,Katzner,S.,Saleem,A. B.,SchoLinck,M. L.,
%et al. (2011). The Detection of Visual Contrast in the Behaving Mouse. Journal
%of Neuroscience,31(31),11351?11361. doi:10.1523/JNEUROSCI.6689-10.2011

%sl=1./(R0 + Rmax.*((c^n)./(C50 + (c^n))))
%This model is more general than the previous model.

%Assign the free parameters to each condition
for i=1:numel(udata)
    
    %Model width of the likelihood
    sl(i)=1./(R0 + Rmax.*((c(i).^n)./(C50 + (c(i).^n))));
    
    %Set prior conditions
    if spe(i)==80
        sp(i)=sp1;
    end
    
    if spe(i)==40
        sp(i)=sp2;
    end
    
    if spe(i)==20
        sp(i)=sp3;
    end
    
    if spe(i)==10
        sp(i)=sp4;
    end
    
    %Bayesian inference
    %mean of the posterior
    uPo(i)=uLl(i).*(1./(1 + (sl(i)./sp(i)).^2)) + (1./(1 + (sp(i)./sl(i)).^2)).*uPr;
    %std of the posterior
    sPo(i)=sqrt(1./((1./sl(i).^2) + (1./sp(i).^2)));
    %std of the estimate
    sEs(i)=sPo(i) + sM;
    
end

%Calculate the SSE
SSE=nansum( ([udata;sdata] - [uPo';sEs']).^2 );

%fprintf('%6.2f  %12.8f\n',[udata;sdata])
%fprintf('%6.2f  %12.8f\n',[uPo;sPo]);
%
%Store the model's informations
pred.mean=uPo;
pred.std=sEs;
fitP=[R0;Rmax;n;C50;sp1;sp2;sp3;sp4;sM];
function [SSE,pred,fitP]=makeSSE3(udata,Fbp,fitPra)
%.........................................................................................................................
%Representations
%.........................................................................................................................-
%c     :coherence (e.g.,12,35,100%).
%uLl   :displayed direction.
%uPr   :prior's mean,i.e.,most likely direction for gaussians.
%uPo   :posterior's mean,i.e.,estimated direction.
%sl    :std of the likelihood.
%sPr   `: std of the prior (e.g.,20 and inf).
%s     :a constant.
%k:    :a free parameter (k=s/sPr)

%Task parameters
%Coherence
c =Fbp(:,1);%coherence
sp=Fbp(:,2);%std of the prior
%Mean of the prior
uPr=225;
%Mean of the likelihood
uLl=Fbp(:,3);
%Free parameter
sl=fitPra;

%Run Bayesian inference to predict the mean of the data
%https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
uPo=(1./( 1+(sl./sp).^2 )).*uLl + (1./( 1+(sp./sl).^2 )).*uPr;
%sPo=1./((1./sl.^2) + (1./sp.^2));

%Calculate the SSE
SSE=nansum( (udata - uPo).^2 );

%Store the model's informations
pred.mean=uPo;
fitP=sl;
data    =udata;
function [SSE,pred,fitP]=makeSSE4(udata,sdata,Fbp,fitPra)

%Representations
%llh
%mean
uLl=Fbp(:,3);
%std
sl=fitPra;

%prior
%mean
uPr=225;
%std
sp=Fbp(:,2);

%Bayes infer
%https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
%mean
uPo=(1./(1+(sl./sp).^2)).*uLl+(1./(1+(sp./sl).^2)).*uPr;
%std
sPo=sqrt(1./((1./sl.^2)+(1./sp.^2)));

%Calculate the SSE,predictions and fit parameters
SSE=nansum(([udata;sdata]-[uPo;sPo]).^2);
pred.mean=uPo;
pred.std=sPo;
fitP=sl;
function [logL,fitP]=makeSSE5(data,displ,coh,Pstd,fitP)
%This nethod is 0.20 seconds faster per function evaluation than 1)
%re-calculating the solution for each data point in a loop and 2) >20
%seconds faster per function evaluation than using fsolve to find the
%solution for each data point in the loop.

%This code can be improved look at GirshickML code in Modeling folder. We
%can get rid of the rounding of the data and of the for loop!


%fit parameters(k has no unit, ideally k is in the range 0:inf)
kl1=fitP(1);
kl2=fitP(2);
kl3=fitP(3);
kp1=fitP(4);
kp2=fitP(5);
kp3=fitP(6);
kp4=fitP(7);
Pm=fitP(8);

%degrees
ulall=1:1:360;
up=225;
ul=nan(numel(data),1);

%3coh*4priors Girshick lookup tables (degrees)
ul11=lookuptable(ulall,up,kl1,kp1);
ul12=lookuptable(ulall,up,kl1,kp2);
ul13=lookuptable(ulall,up,kl1,kp3);
ul14=lookuptable(ulall,up,kl1,kp4);
ul21=lookuptable(ulall,up,kl2,kp1);
ul22=lookuptable(ulall,up,kl2,kp2);
ul23=lookuptable(ulall,up,kl2,kp3);
ul24=lookuptable(ulall,up,kl2,kp4);
ul31=lookuptable(ulall,up,kl3,kp1);
ul32=lookuptable(ulall,up,kl3,kp2);
ul33=lookuptable(ulall,up,kl3,kp3);
ul34=lookuptable(ulall,up,kl3,kp4);

%assign ul solution (degrees) to each data point
data=round(data);
%make sure data ranges in 1:1:360.I saw that simulateData produces both 0
%and 360 that are the same. All 0 will be 360.
%parpool parfor if possible
data(data==0)=360;
for i=1:360
    ul(data==i&coh==0.24&Pstd==80)=ul11(i);
    ul(data==i&coh==0.24&Pstd==40)=ul12(i);
    ul(data==i&coh==0.24&Pstd==20)=ul13(i);
    ul(data==i&coh==0.24&Pstd==10)=ul14(i);
    ul(data==i&coh==0.12&Pstd==80)=ul21(i);
    ul(data==i&coh==0.12&Pstd==40)=ul22(i);
    ul(data==i&coh==0.12&Pstd==20)=ul23(i);
    ul(data==i&coh==0.12&Pstd==10)=ul24(i);
    ul(data==i&coh==0.06&Pstd==80)=ul31(i);
    ul(data==i&coh==0.06&Pstd==40)=ul32(i);
    ul(data==i&coh==0.06&Pstd==20)=ul33(i);
    ul(data==i&coh==0.06&Pstd==10)=ul34(i);
end

%The likelihood of observing data i is the probability that measurement
%mi, that is the mean of the likelihood (ul(i)), is observed given a
%measurement distribution which mean is the displayed direction + the
%likelihood of observing this measurement if it comes from motor noise "Pm".
m=1:1:360;
um=displ;
kl(coh==0.24)=fitP(1);
kl(coh==0.12)=fitP(2);
kl(coh==0.06)=fitP(3);
mpdf=vmPdfs(m,um,kl);

%likelihood
mpdf

Pm=0;% to remove
Pe=(1-Pm)*mpdf+Pm;



%single trial's measurement, its position(row) for each trial(col) and its
%probability (also likelihood of trial's data).
mi=ul;
mipos=mi+numel(m)*((0:1:numel(mi)-1)');
ML.p=mpdf(mipos);

%We use log likelihood because likelihood is so small that matlab cannot
%encode it properly (numerical unstability).
logL=-sum(log(ML.p));

%Look at fitting
fprintf('\n %4.2f %4.2g %4.2g %4.2g %4.2g %4.2g %4.2g %4.2g %4.2g \n',logL,kl1,kl2,...
    kl3,kp1,kp2,kp3,kp4,Pm)
pred=makePre5(displ,coh,Pstd,[kl1 kl2 kl3 kp1 kp2 kp3 kp4]);
drawPred5(pred,data,displ,coh,Pstd);
drawnow
%logl
% subplot(3,1,1)
% title(sprintf('%s %f','logL:',logL));
% hold all;
% hi=plot(1,1,'.w');
% plot(hi,logL,'o','Markerfacecolor','k')
% drawnow
% %kl
% subplot(3,1,2)
% title(sprintf('%s %6.4f %s %6.4f %s %6.4f',...
%     ' kl1: ',kl1,' kl2: ',kl2,' kl3: ',kl3));
% hold all;
% plot(hi,kl1,'o','Markerfacecolor','r')
% plot(hi,kl2,'o','Markerfacecolor','b')
% plot(hi,kl3,'o','Markerfacecolor','k')
% drawnow
% %kp
% subplot(3,1,3)
% title(sprintf('%s %6.4f %s %6.4f %s %6.4f %s %6.4f',...
%     ' kp1: ',kp1,' kp2: ',kp2,' kp3: ',kp3,' kp4: ',kp4));
% hold all;
% plot(hi,kp1,'o','Markerfacecolor','r')
% plot(hi,kp2,'o','Markerfacecolor','b')
% plot(hi,kp3,'o','Markerfacecolor','k')
% plot(hi,kp4,'o','Markerfacecolor','k')
% drawnow
function [SSE,pred,fitP]=makeSSE6(udata,sdata,Fbp,fitP)

%Llh
%mean&std
uLl=Fbp(:,3);
sl=nan(numel(udata),1);
sl(Fbp(:,1)==0.24)=fitP(1);
sl(Fbp(:,1)==0.12)=fitP(2);
sl(Fbp(:,1)==0.06)=fitP(3);


%Prior
%mean&std
uPr=225;
sp=nan(numel(udata),1);
sp(Fbp(:,2)==80)=fitP(4);
sp(Fbp(:,2)==40)=fitP(5);
sp(Fbp(:,2)==20)=fitP(6);
sp(Fbp(:,2)==10)=fitP(7);

%Motor noise
sM=fitP(8);


%Posterior
%mean&std
%uPo=(1./(1+(sl./sp).^2)).*uLl+(1./(1+(sp./sl).^2)).*uPr;
%sPo=sqrt(1./((1./sl.^2)+(1./sp.^2)));

%Posterior
%Ideally calculating mean and std requires integrating over [-inf:+inf].
%Here we only get approximations.
%!!! the cost of this loop can be reduced by looping over unique(uLl).
%note: when sl~0 the equation stay stuck at inf*0.
x=[-1000:10:1000]';
for i=1:numel(uLl)
    %llh
    l(:,i)=((sl(i)*sqrt(2*pi)).^-1).*exp(-0.5*((x-uLl(i)).^2)/(sl(i).^2));
    l(:,i)=l(:,i)/sum(l(:,i));
    %note: approximation of uL is very bad but gaussian is ok.
    
    %prior
    pr(:,i)=((sp(i)*sqrt(2*pi)).^-1).*exp(-0.5*((x-uPr).^2)/(sp(i).^2));
    pr(:,i)=pr(:,i)/sum(pr(:,i));
    
    %posterior
    po(:,i)=l(:,i).*pr(:,i)/sum(l(:,i).*pr(:,i));
    uPo(i,1)=sum(x.*po(:,i));
    sPo(i,1)=sum(sqrt(po(:,i).*(x-uPo(i)).^2));
    
    %estimate
    sEs(i,1)=sPo(i)+sM;
end


%SSE,predictions and parameters
SSE=nansum(([udata;sdata]-[uPo;sEs]).^2);
pred.mean=uPo;
pred.std=sEs;


%R^2
function R2=makeR2(data,SSE)
%R2=1 - SSE/SST
%SSE (alias residual) is the SSE calculated for the best model's prediction.
%note: R^ (alias coefficient of determination) can be <0 when the model
%does an awful job predicting the data. In this case SSE exceeds that is
%the model fits the data even worse than does a horizontal line.
%Model may not be appropriate or constraints may not be set correctly.
%see http://www.graphpad.com/support/faqid/711/
%
SST=sum((data - mean(data)).^2);
R2=1 - SSE/SST;
function R2=makeR22(data,pred)
%R2=1 - SSE/SST
%SSE (alias residual) is the SSE calculated for the best model's prediction.
%note: R^ (alias coefficient of determination) can be <0 when the model
%does an awful job predicting the data. In this case SSE exceeds that is
%the model fits the data even worse than does a horizontal line.
%Model may not be appropriate or constraints may not be set correctly.
%see http://www.graphpad.com/support/faqid/711/
%
R=corr(data,pred);
R2=R^2;



%Make predictions
function [pred,F]=makePre1(F,BestfitP)

%fixed parameters
c=F.g1L;
uPr=225;
uLl=F.g3L;

%fit parameters
%k=cte/width prior;
k=BestfitP;

%Operations
%# in theory,s/c is the variance of the likelihood (i.e.,noise) and can't
%be negative
%# in theory sPr is the variance of the prior and can't be negative.
%# What about k=s/sPr ?

%sPr cannot < 0,thus if k<0,it means s<0. If s<0,s/c>0 only if c<0;
%Thus,when k<0,encoding of coherence<0;
%Not sure how that makes sense....

%Bayesian inference
for i=1:numel(F.g1L)
    for j=1:numel(F.g2L)
        for l=1:numel(F.g3L)
            uPo(l,i,j)=(1./(1+(k(j)./c(i)).^2)).*uLl(l)+uPr.*(1./(1+(c(i)./k(j)).^2));
        end
    end
end

%Store the model's informations
pred.mean=uPo;
function [pred,F]=makePre2a(F,BestfitP)
%.........................................................................................................................
%Representations
%.........................................................................................................................-
%c     :coherence (e.g.,12,35,100%).
%uLl   :displayed direction.
%uPr   :prior's mean,i.e.,most likely direction for gaussians.
%uPo   :posterior's mean,i.e.,estimated direction.
%sl    :std of the likelihood.
%sPr   `: std of the prior (e.g.,20 and inf).
%s     :a constant.

%udata :
%sdata :

%Representations
%.........................................................................................................................-
%Task parameters
c=F.g1L;%coherence
uPr=225;%Mean of the prior
uLl=F.g3L;%Mean of the likelihood

%Free parameters
%llh
R0=BestfitP(1);
Rmax=BestfitP(2);
n=BestfitP(3);
C50=BestfitP(4);
%prior
sp(1)=BestfitP(5);%std of prior 1
sp(2)=BestfitP(6);
sp(3)=BestfitP(7);
sp(4)=BestfitP(8);
%motor
sM=BestfitP(9);%motor noise

%Debuging
%fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
%    [sl1;sl2;sl3;sp1;sp2;sp3;sp4])

%Assign the free parameters to each condition
%loop over
for i=1:numel(F.g1L)%e.g.,coherence
    for j=1:numel(F.g2L) %e.g.,prior
        for k=1:numel(F.g3L)
            
            %Run Bayesian inference
            %https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
            
            %width of the likelihood
            sl(i)=1./(R0 + Rmax.*((c(i).^n)./(C50 + (c(i).^n))));
            
            %mean of the posterior
            uPo(k,i,j)=(1./( 1+(sl(i)./sp(j)).^2 )).*uLl(k) + (1./( 1+(sp(j)./sl(i)).^2 )).*uPr;
            %std of the posterior
            sPo(k,i,j)=sqrt(1./((1./sl(i).^2) + (1./sp(j).^2)));
            %std of the estimate
            sEs(k,i,j)=sPo(k,i,j) + sM;
            
            %Store conditions
            F.g1bp(k,i,j)=F.g1L(i);
            F.g2bp(k,i,j)=F.g2L(j);
            F.g3bp(k,i,j)=F.g3L(k);
        end
    end
end

%Store the model's informations
pred.mean=uPo;
pred.std=sEs;
function [pred,F]=makePre(F,BestfitP)
%.........................................................................................................................
%Representations
%.........................................................................................................................-
%c     :coherence (e.g.,12,35,100%).
%uLl   :displayed direction.
%uPr   :prior's mean,i.e.,most likely direction for gaussians.
%uPo   :posterior's mean,i.e.,estimated direction.
%sl    :std of the likelihood.
%sPr   `: std of the prior (e.g.,20 and inf).
%s     :a constant.

%udata :
%sdata :

%Task parameters
uPr=225;%mean of the prior

%Input
uLl=F.g3L;

%Set the free parameters
sl(1)=BestfitP(1);%std of likelihood 1
sl(2)=BestfitP(2);
sl(3)=BestfitP(3);
sp(1)=BestfitP(4);%std of prior 1
sp(2)=BestfitP(5);
sp(3)=BestfitP(6);
sp(4)=BestfitP(7);
sM =BestfitP(8);%motor noise

%Debuging
%fprintf('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
%    [sl1;sl2;sl3;sp1;sp2;sp3;sp4])

%Assign the free parameters to each condition
%loop over
for i=1:numel(F.g1L)%e.g.,coherence
    for j=1:numel(F.g2L) %e.g.,prior
        for k=1:numel(F.g3L)
            
            %Run Bayesian inference
            %https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
            %mean of the posterior
            uPo(k,i,j)=(1./( 1+(sl(i)./sp(j)).^2 )).*uLl(k) + (1./( 1+(sp(j)./sl(i)).^2 )).*uPr;
            %std of the posterior
            sPo(k,i,j)=sqrt(1./((1./sl(i).^2) + (1./sp(j).^2)));
            %std of the estimate
            sEs(k,i,j)=sPo(k,i,j) + sM;
            
            
            %Store conditions
            F.g1bp(k,i,j)=F.g1L(i);
            F.g2bp(k,i,j)=F.g2L(j);
            F.g3bp(k,i,j)=F.g3L(k);
        end
    end
end

%Store the model's informations
pred.mean=uPo;
pred.std=sEs;
function [pred]=makePre5(d,coh,pstd,k)
% tic
%Create a measurement density at each trial which mean is the
%displayed direction and width is km=kl, the fit parameter calculated for
%this trial.

%Sample the measurement density and get a measurement mi for this
%trial. For trial-average prediction mi="displayed direction";

%mi=ul, likelihood's mean in this trial. Input ul, up(known), kl
%and kp (best fit parameters) in the equation that calculate the mean
%of the posterior (Murray and Morgenstern) for this trial, that is the
%trial predicted estimate.

%measurement densities and priors' strengths
up=225;
kml(coh==0.24)=k(1);
kml(coh==0.12)=k(2);
kml(coh==0.06)=k(3);
kp(pstd==80)=k(4);
kp(pstd==40)=k(5);
kp(pstd==20)=k(6);
kp(pstd==10)=k(7);

%estimates
xe=1:1:360;

%%trial-average prediction
%m=vmPdfs(xe,d,kml);
mi=d;
pred=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)),(kml./kp)'+cos(de2r(up,1)-de2r(mi,1))))';

%measurement densities
% m=vmPdfs(xe,d,kml);
% mi=nan(1,numel(d));
% for i=1:numel(d)
%
%    %sample
%    mi(i)=randsample(xe,1,true,m(:,i));
% end
% %trial-estimate.
% pred=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)), (kml./kp)+cos(de2r(up,1)-de2r(mi,1))));

% toc

%make sure data ranges in 1:1:360.
pred(pred==0)=360;


%Draw predictions
function [fig1]=drawPred(pred,udata,Fbp,F,fig)
global databank

%INPUT:
%F: conditions
%udata: mean data
%pred: predictions

%Draw group1-level 1 (e.g.,std prior 80)
%Enlarge figure for good quality publication
%fig1.hdle=figure('Position',[0 0 1000 400]);%pixels
fig1.hdle=figure('color',[1 1 1]);
fig1.nm=[fig.nm,'_DataAndModel'];

%Set factors
%F 1 (e.g.,coherence)
clear i
for i=1:F.g1.numL
    indX.g1.Ll_i(i)={find(Fbp(:,1)==F.g1.L(i))};
end
%F 2 (e.g.,priors)
clear i
for i=1:F.g2.numL
    indX.g2.Ll_i(i)={find(Fbp(:,2)==F.g2.L(i))};
end
%F 3 (e.g.,displayed directions)
clear i
for i=1:F.g3.numL
    indX.g3.Ll_i(i)={find(Fbp(:,3)==F.g3.L(i))};
end

%initialize the parameters of the plot
%set the colors
c=colormap;
F.g2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0],...
    c(40,:),...
    c(32,:),...model
    c(27,:),...
    c(22,:),...
    c(8 ,:),...
    c(5 ,:),...
    c(1 ,:)};%group 2(e.g.,coherences)

if numel(F.g2.color) < F.g2.numL
    disp(['--- You may want to add ',...
        num2str(F.g2.numL - numel(F.g2.color)),...
        'more colors to the color code ---']);
    return
end
%set axes positions
width=1/(F.g1.numL+1);
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end
%.................
%Plot data - mean
%.................
%space axes
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

%loop over levels of factor 1 (e.g. coherences)
for k=1:F.g1.numL
    ax(k)=axes('position',axs.position(k,:));
    axis square
    
    %priors' mean
    hold all
    p1=plot([F.g3.L(1) F.g3.L(end)],[databank.data{1,7} databank.data{1,7}],...
        'b:',...
        'linewidth',1,...
        'DisplayName','Prior mean');
    
    %ideal performance
    p10=plot(F.g3.L,F.g3.L','k:',...
        'linewidth',1,...
        'Displayname','Ideal predictions');
    
    for j=1:F.g2.numL
        %positions of conditions
        indX.g1g2(j,k)={intersect( indX.g1.Ll_i{k},indX.g2.Ll_i{j} ) };
        %store data
        fig1.datamean{j,k}=udata(indX.g1g2{j,k});%es{i,j,k}.deg.mean;
        fig1.dataInfog1{j,k}=Fbp(indX.g1g2{j,k},1);%F.g1.L(k);
        fig1.dataInfog2{j,k}=Fbp(indX.g1g2{j,k},2);%F.g2.L(j);
        fig1.dataInfog3{j,k}=Fbp(indX.g1g2{j,k},3);%F.g3.L(i);
        
        
        %        %Normalize data for plotting (because circular data) -already
        %        done earlier in the code. should be suppressed
        %        %..............................................................................................................----
        %        %In case the direction is displayed in the 3rd quarter.
        %        if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
        %            %"Linearize" the value of the estimated direction
        %            fig1.datamean{j,k}=360 + fig1.datamean{j,k};
        %            %In case the direction is displayed in the 1st quarter.
        %        elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
        %            %"Linearize" the value of the estimated direction
        %            fig1.datamean{j,k}=fig1.datamean{j,k} - 360;
        %        end
        
        
        %data
        scatter(fig1.dataInfog3{j,k},fig1.datamean{j,k},...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.g2.color{j},...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))))
        %p12=[p12 p12_];
        
        %predictions
        p13a=plot(F.g3.L,pred.mean(:,k,j),'-',...
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',2,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
    end
    
    %graphics
    %x-unit
    xunit=1:6:F.g3.numL;
    set(gca,...
        'xtick',F.g3.L(xunit),'xticklabel',F.g3.L(xunit),...
        'ytick',F.g3.L(xunit),'yticklabel',F.g3.L(xunit),...
        'fontsize',12);
    xlim([min(F.g3.L)-10 max(F.g3.L)+10]);
    ylim([min(F.g3.L)-10 max(F.g3.L)+10]);
    if k==1
        ylabel(ax(k),'Estimated directions','fontsize',12);
    end
    xlabel('Displayed directions','fontsize',12);
    drawPublishAxis
end
%%Legend the graph
%leg=legend (p12,'location','BestOutside');
%set(leg,'Box','off');

%Title
ax(k+1,:)=axes('position',[0 0 0.97 0.75],'visible','off');
for k=1:F.g1.numL
    %    text(axs.position(k,1)+axs.position(k,3)/2,...
    %        0.05 + axs.position(k,2)+axs.position(k,4),...
    %        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
    %        'fontweight','Bold',...
    %        'fontsize',12);
    xpos(k)=axs.position(k,1) + axs.position(k,3)/2;
    ypos(k)=0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
        'fontweight','Bold',...
        'fontsize',12);
    
end
function [fig1]=drawPred4(pred,udata,sdata,Fbp,F,fig)
global databank

%Draw group1-Ll1 (e.g.,prior std=80)
%Enlarge figure for good quality publication
%fig1.hdle=figure('Position',[0 0 1000 400]);%pixels
fig1.hdle=figure('color',[1 1 1]);
fig1.nm=[fig.nm,'_DataAndModel'];

%Organize the data for plotting
Fbp4plot=[];
udata4plot   =[];
sdata4plot   =[];
pred4plotmean=[];
pred4plotstd=[];

for i=1:numel(Fbp)
    %mean of the data (e.g.,estimated directions)
    udata4plot=[udata4plot;udata{i}];
    %sdt of the data (e.g.,std estimated directions)
    sdata4plot=[sdata4plot;sdata{i}];
    %model's predictions mean and std
    pred4plotmean=[pred4plotmean;pred{i}.mean];
    pred4plotstd=[pred4plotstd;pred{i}.std];
    %task parameters
    Fbp4plot=[Fbp4plot;Fbp{i}];
end

%%Set factors
%F 1 (e.g.,coherence)
clear i
for i=1:F.g1.numL
    indX.g1.Ll_i(i)={find(Fbp4plot(:,1)==F.g1.L(i))};
end
%F 2 (e.g.,priors)
clear i
for i=1:F.g2.numL
    indX.g2.Ll_i(i)={find(Fbp4plot(:,2)==F.g2.L(i))};
end
%F 3 (e.g.,displayed directions)
clear i
for i=1:F.g3.numL
    indX.g3.Ll_i(i)={find(Fbp4plot(:,3)==F.g3.L(i))};
end


%initialize the parameters of the plot
%set the colors
c=colormap;
F.g2.color={[0.5 0 0],...
    c(55,:),...
    [1 0.4 0],...
    [0.75 0.75 0],...
    c(40,:),...
    c(32,:),...
    c(27,:),...
    c(22,:),...
    c(8 ,:),...
    c(5 ,:),...
    c(1 ,:)};%group 2(e.g.,coherences)

if numel(F.g2.color) < F.g2.numL
    disp(['--- You may want to add ',...
        num2str(F.g2.numL - numel(F.g2.color)),...
        'more colors to the color code ---']);
    return
end
%set axes positions
width=1/(F.g1.numL+1);
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end
%.........................................................................................................................----
%Plot data - mean
%.........................................................................................................................----
%space axes
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

%loop over levels of factor 1 (e.g. coherences)
for k=1:F.g1.numL
    ax(k)=axes('position',axs.position(k,:));
    axis square
    
    %draw background
    %.........................................................................................................................---
    hold all
    %draw priors' mean
    p1=plot([F.g3.L(1) F.g3.L(end)],[databank.data{1,7} databank.data{1,7}],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    %Note: add linecode here for drawing an arrow indicating initial position
    %of response line.
    
    %draw the ideal performance line
    p10=plot([F.g3.L],[F.g3.L'],'k:',...
        'linewidth',1.00005,...
        'Displayname','Ideal predictions');
    
    %organize data for plotting
    %.........................................................................................................................---
    %change the level of factor 2 (e.g.,prior's strength)
    for j=1:F.g2.numL
        %change the level of factor 3 (displayed directions)
        %coordinates of single conditions
        indX.g1g2(j,k)={intersect( indX.g1.Ll_i{k},indX.g2.Ll_i{j} ) };
        %store data
        fig1.datamean{j,k}=udata4plot(indX.g1g2{j,k});%es{i,j,k}.deg.mean;
        
        %extract information about data
        %g1
        fig1.dataInfog1{j,k}=Fbp4plot(indX.g1g2{j,k},1);%F.g1.L(k);
        %g2
        fig1.dataInfog2{j,k}=Fbp4plot(indX.g1g2{j,k},2);%F.g2.L(j);
        %g3
        fig1.dataInfog3{j,k}=Fbp4plot(indX.g1g2{j,k},3);%F.g3.L(i);
        
        %store model's predictions
        fig1.pred.mean{j,k}=pred4plotmean(indX.g1g2{j,k});
        fig1.pred.std{j,k}=pred4plotstd(indX.g1g2{j,k});
        
        
        %        %Normalize data for plotting (because circular data) -already
        %        done earlier in the code. should be suppressed
        %        %..............................................................................................................----
        %        %In case the direction is displayed in the 3rd quarter.
        %        if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
        %            %"Linearize" the value of the estimated direction
        %            fig1.datamean{j,k}=360 + fig1.datamean{j,k};
        %            %In case the direction is displayed in the 1st quarter.
        %        elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
        %            %"Linearize" the value of the estimated direction
        %            fig1.datamean{j,k}=fig1.datamean{j,k} - 360;
        %        end
        
        
        %plot data
        %..............................................................................................................-----
        p12_=plot(fig1.dataInfog3{j,k},fig1.datamean{j,k},'o',... %groups 3 forms x-axis
            'color',F.g2.color{j},...
            'markerfacecolor',F.g2.color{j},...
            'markersize',8,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
        %p12=[p12 p12_];
        
        
        %plot model's predictions
        %..............................................................................................................----
        %mean
        p13a=plot(fig1.dataInfog3{j,k} ,[fig1.pred.mean{j,k}],'-',... %groups 3 forms x-axis
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',1.0005,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
    end
    
    %set the graph parameters
    %set the unit steps of the x axis
    xunit=1:6:F.g3.numL;
    set(gca,...
        'xtick',F.g3.L(xunit),'xticklabel',F.g3.L(xunit),...
        'ytick',F.g3.L(xunit),'yticklabel',F.g3.L(xunit),...
        'fontsize',12);
    %Set the limits of the axes
    xlim([min(F.g3.L)-10 max(F.g3.L)+10]);
    ylim([min(F.g3.L)-10 max(F.g3.L)+10]);
    %set the labels of the axes
    if k==1
        ylabel(ax(k),'Estimated directions','fontsize',12);
    end
    xlabel('Displayed directions','fontsize',12);
end
%%Legend the graph
%leg=legend (p12,'location','BestOutside');
%set(leg,'Box','off');

%Title
ax(k+1,:)=axes('position',[0 0 0.97 0.75],'visible','off');
for k=1:F.g1.numL
    %    text(axs.position(k,1)+axs.position(k,3)/2,...
    %        0.05 + axs.position(k,2)+axs.position(k,4),...
    %        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
    %        'fontweight','Bold',...
    %        'fontsize',12);
    xpos(k)=axs.position(k,1) + axs.position(k,3)/2;
    ypos(k)=0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
        'fontweight','Bold',...
        'fontsize',12);
    
end


drawPublishAxis

%.........................................................................................................................----
%Plot data - std
%.........................................................................................................................----

fig2.hdle=figure('color','w');
fig2.nm=[fig.nm,'_std'];

%set the axes positions
width=1/(F.g1.numL+1);

%space axes
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

%loop over levels of factor 1 (e.g. coherences)
for k=1:F.g1.numL
    ax(k)=axes('position',axs.position(k,:));
    axis square
    
    %draw background
    %.........................................................................................................................---
    hold all
    %draw priors' mean
    p1=plot([databank.data{1,7} databank.data{1,7}],[0 150],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    %organize data for plotting
    %.........................................................................................................................---
    %change the level of factor 2 (e.g.,prior's strength)
    for j=1:F.g2.numL
        %change the level of factor 3 (displayed directions)
        %coordinates of single conditions
        indX.g1g2(j,k)={intersect( indX.g1.Ll_i{k},indX.g2.Ll_i{j} ) };
        
        %store data
        fig1.datastd{j,k}=sdata4plot(indX.g1g2{j,k});
        
        %extract information about data
        %g1
        fig1.dataInfog1{j,k}=Fbp4plot(indX.g1g2{j,k},1);%F.g1.L(k);
        %g2
        fig1.dataInfog2{j,k}=Fbp4plot(indX.g1g2{j,k},2);%F.g2.L(j);
        %g3
        fig1.dataInfog3{j,k}=Fbp4plot(indX.g1g2{j,k},3);%F.g3.L(i);
        
        %Store model's predictions
        fig1.pred.std{j,k}=pred4plotstd(indX.g1g2{j,k});
        
        
        %        %Normalize data for plotting (because circular data) -already
        %        done earlier in the code. should be suppressed
        %        %..............................................................................................................----
        %        %In case the direction is displayed in the 3rd quarter.
        %        if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
        %            %"Linearize" the value of the estimated direction
        %            fig1.datamean{j,k}=360 + fig1.datamean{j,k};
        %            %In case the direction is displayed in the 1st quarter.
        %        elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
        %            %"Linearize" the value of the estimated direction
        %            fig1.datamean{j,k}=fig1.datamean{j,k} - 360;
        %        end
        
        
        %plot data
        %..............................................................................................................-----
        p12_=plot(fig1.dataInfog3{j,k},fig1.datastd{j,k},'o',... %groups 3 forms x-axis
            'color',F.g2.color{j},...
            'markerfacecolor',F.g2.color{j},...
            'markersize',8,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
        
        %plot model's predictions
        %..............................................................................................................----
        %mean
        p13a=plot(fig1.dataInfog3{j,k} ,[fig1.pred.std{j,k}],'-',... %groups 3 forms x-axis
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',1.0005,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
    end
    
    %set the graph parameters
    %set the unit steps of the x axis
    xunit=1:6:F.g3.numL;
    set(gca,...
        'xtick',F.g3.L(xunit),'xticklabel',F.g3.L(xunit),...
        'fontsize',12);
    %Set the limits of the axes
    xlim([min(F.g3.L)-10 max(F.g3.L)+10]);
    ylim([0 150]);
    %set the labels of the axes
    if k==1
        ylabel(ax(k),'Std of estimated directions','fontsize',12);
    end
    xlabel('Displayed directions','fontsize',12);
end
%%Legend the graph
%leg=legend (p12,'location','BestOutside');
%set(leg,'Box','off');

%Title
ax(k+1,:)=axes('position',[0 0 0.97 0.75],'visible','off');
for k=1:F.g1.numL
    %    text(axs.position(k,1)+axs.position(k,3)/2,...
    %        0.05 + axs.position(k,2)+axs.position(k,4),...
    %        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
    %        'fontweight','Bold',...
    %        'fontsize',12);
    xpos(k)=axs.position(k,1) + axs.position(k,3)/2;
    ypos(k)=0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
        'fontweight','Bold',...
        'fontsize',12);
    
end


%backup figures informations
figs={fig1,fig2};
function [meansD,meansP,stdsD,stdsP]=drawPred5(pred,data,d,coh,pstd)

%data(cartesians)
datac=polar2cartesian(data,1);
predc=polar2cartesian(pred',1);

%factors 1,2,3
F.f1.i=d;
F.f1.nm='d';
F.f1.L=unique(F.f1.i);
F.f1.L=sort(F.f1.L,'ascend');
F.f1.n=numel(F.f1.L);

F.f2.i=coh;
F.f2.nm='coh';
F.f2.L=unique(F.f2.i);
F.f2.L=sort(F.f2.L,'descend');
F.f2.n=numel(F.f2.L);

F.f3.i=pstd;
F.f3.nm='pstd';
F.f3.L=unique(F.f3.i);
F.f3.L=sort(F.f3.L,'descend');
F.f3.n=numel(F.f3.L);

%positions main
for i=1:F.f1.n
    F.f1.pos(i)={find(F.f1.i==F.f1.L(i))};
end
for i=1:F.f2.n
    F.f2.pos(i)={find(F.f2.i==F.f2.L(i))};
end
for i=1:F.f3.n
    F.f3.pos(i)={find(F.f3.i==F.f3.L(i))};
end

%positions inter
for k=1:F.f1.n
    for j=1:F.f2.n
        for i=1:F.f3.n
            F.inter.pos(k,i,j)=...
                {intersect( ...
                intersect(F.f1.pos{k},F.f2.pos{j}),...
                F.f3.pos{i})};
        end
    end
end

%Make mean & std
c=colormap;
F.f2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};
for j=1:F.f2.n
    for k=1:F.f1.n
        for i=1:F.f3.n
            %data
            D{k,i,j}=statcircular(datac(F.inter.pos{k,i,j},:));
            meansD(k,i,j)=D{k,i,j}.deg.mean;
            
            %pred
            P{k,i,j}=statcircular(predc(F.inter.pos{k,i,j},:));
            meansP(k,i,j)=P{k,i,j}.deg.mean;
            
            stdsD(k,i,j)=D{k,i,j}.deg.std;
            stdsP(k,i,j)=P{k,i,j}.deg.std;
        end
    end
end

%mean
h(1)=figure(1);
clf
set(h(1),'color','w')
ti={'Coh: 24%','Coh: 12%','Coh: 6%'};
for j=1:F.f2.n
    subplot(1,3,j)
    for i=1:F.f3.n
        hold all
        %data
        scatter(F.f1.L(:),meansD(:,i,j)',...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.f2.color{i},...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        %pred
        plot(F.f1.L(:),meansP(:,i,j)','-',...
            'color',F.f2.color{i},...
            'linewidth',2,...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        %ideal
        plot(F.f1.L(:),F.f1.L(:),'k.','Markersize',3)
    end
    xlim([0 360])
    ylim([0 360])
    xlabel('Displayed direction (?)')
    ylabel('Estimated direction (?)')
    title(ti{j})
    drawPublishAxis
end

% %std
% h(2)=figure(2);
% clf
% set(h(2),'color','w')
% for j=1:F.f2.n
%     subplot(1,3,j)
%     for i=1:F.f3.n
%         hold all
%         %data
%         scatter(F.f1.L(:),stdsD(:,i,j)',...
%             'MarkerEdgeColor','w',...
%             'MarkerFaceColor',F.f2.color{i},...
%             'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
%         %pred
%         plot(F.f1.L(:),stdsP(:,i,j)','-',...
%             'color',F.f2.color{i},...
%             'linewidth',2,...
%             'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
%     end
%     xlim([0 360])
%     ylim([0 360])
%     xlabel('Displayed direction (?)')
%     ylabel('Estimated direction (?)')
%     title(ti{j})
%     drawPublishAxis
% end
function [fig1]=drawPred6(pred,F,fig,datafitt)
global databank

%set figure
fig1.hdle=figure('color',[1 1 1]);
fig1.nm=[fig.nm,'_DataAndModel'];

%set data and factors to plot
udata=datafitt.data{strcmp(datafitt.nm,'mean')==1};
sdata=datafitt.data{strcmp(datafitt.nm,'std')==1};
F.g1.Fbp=datafitt.data{strcmp(datafitt.nm,'F1')==1};
F.g2.Fbp=datafitt.data{strcmp(datafitt.nm,'F2')==1};
F.g3.Fbp=datafitt.data{strcmp(datafitt.nm,'F3')==1};

%normalized or not
if ~isfield(F.g3,'distL'); F.g3.disp=F.g3.L;
else F.g3.disp=F.g3.distL; end

%Graphics
%set colors
c=colormap;
F.g2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0],...
    c(40,:),...
    c(32,:),...
    c(27,:),...
    c(22,:),...
    c(8 ,:),...
    c(5 ,:),...
    c(1 ,:)};
%warning
if numel(F.g2.color)<F.g2.numL
    disp(['--- You may want to add ',...
        num2str(F.g2.numL - numel(F.g2.color)),...
        'more colors to the color code ---']);
    return
end

%make axes
axs=makeAxes(F);

%Draw mean
%factor 1(coh)
for k=1:F.g1.numL
    %axe
    ax(k)=axes('position',axs.position(k,:));
    axis square
    
    %priors' mean
    hold all
    p1=plot([F.g3.disp(1) F.g3.disp(end)],[225 225],...
        'k:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    %ideal predictions
    p10=plot(F.g3.disp,F.g3.disp','k:',...
        'linewidth',1,...
        'Displayname','Ideal predictions');
    
    %factor 2(prior)
    for j=1:F.g2.numL
        %mean
        scatter(F.g3.disp,udata(:,j,k)',...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.g2.color{j},...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))))
        
        %predictions
        plot(F.g3.disp,pred.mean(:,k,j),'-',... %groups 3 forms x-axis
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',2,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
    end
    
    %set graph
    %set x-unit
    %xunit=1:3:F.g3.numL;
    %set(gca,'xtick',F.g3.disp(xunit),'xticklabel',F.g3.disp(xunit),...
    %   'ytick',F.g3.disp(xunit),'yticklabel',F.g3.disp(xunit),...
    %  'fontsize',12);
    set(gca,'xtick',[0 90 180 270],'xticklabel',[0 90 180 270],...
        'ytick',[0 90 180 270],'yticklabel',[0 90 180 270],...
        'fontsize',12);
    %set axes' limits
    xlim([min(F.g3.disp)-11 max(abs(F.g3.disp))+11]);
    ylim([min(F.g3.disp)-11 max(abs(F.g3.disp))+11]);
    %set axis' labels
    if k==1
        ylabel(ax(k),'Estimated directions','fontsize',12);
    end
    xlabel('Displayed directions','fontsize',12);
end

%Title
ax(k+1,:)=axes('position',[0 0 0.97 0.75],'visible','off');
for k=1:F.g1.numL
    xpos(k)=axs.position(k,1) + axs.position(k,3)/2;
    ypos(k)=0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
        'fontweight','Bold',...
        'fontsize',12);
end





%Draw std
%set figure
fig2.hdle=figure('color','w');
fig2.nm=[fig.nm,'_std'];

%make axes
axs=makeAxes(F);

%factor 1
for k=1:F.g1.numL
    
    %axes
    ax(k)=axes('position',axs.position(k,:));
    axis square
    
    %prior mean
    hold all
    p1=plot([databank.data{1,7} databank.data{1,7}],[0 150],...
        'k:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    for j=1:F.g2.numL
        %std
        scatter(F.g3.disp,sdata(:,j,k)',60,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.g2.color{j},...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))))
        
        %draw predictions
        p13a=plot(F.g3.disp,pred.std(:,k,j)','-',...
            'color',F.g2.color{j}-[0.2 0 0],...
            'linewidth',2,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
    end
    
    %set graph
    %set x units
    xunit=1:6:F.g3.numL;
    %set(gca,...
    %    'xtick',F.g3.disp(xunit),'xticklabel',F.g3.disp(xunit),...
    %    'fontsize',12)
    set(gca,'xtick',[0 90 180 270],'xticklabel',[0 90 180 270],...
        'fontsize',12);
    %Set axes' limits
    xlim([min(F.g3.disp)-11 max(abs(F.g3.disp))+11]);
    ylim([0 150]);
    
    %set axes' labels
    if k==1
        ylabel(ax(k),'Std of estimated directions','fontsize',12);
    end
    xlabel('Displayed directions','fontsize',12);
end

%Title
ax(k+1,:)=axes('position',[0 0 0.97 0.75],'visible','off');
for k=1:F.g1.numL
    xpos(k)=axs.position(k,1) + axs.position(k,3)/2;
    ypos(k)=0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
        'fontweight','Bold',...
        'fontsize',12);
end


%backup figures informations
figs={fig1,fig2};


%std of fit parameters
function stdPa=GetStdP(Hessian,data,pred,fitP)
%note: complexe number Std means parameters have negative variance

%Calculate parameters' covariance matrix
%d.covar=inv(jacobian'*jacobian);
d.covar=inv(Hessian); %a matrix

%Calculate noise variance
residual=(data-pred')';
residual=residual(~isnan(residual));
noiseVariance=(residual*residual')/(numel(data)-numel(fitP));%a scalar

%Calculate std of model parameters
stdPa=diag(sqrt(noiseVariance*d.covar))';


%Internal priors
function Sp=GetPrior1(freeP,stdPa,factor)
%Convert cell input to matrix
if iscell(freeP)==1;freeP=[freeP{:}];end

%Remove cases of infinite prior from the analysis
if freeP(:,1)==+inf;freeP(1)=[];end

%priors' experimental&estimated data
Sp.exp.raw=factor;
Sp.es.raw=freeP;

%std of parameters
if isempty(stdPa)==0;Sp.std=stdPa;else Sp.std=[];end

%Normalized priors' widths
%Experimental&estimated
for i=1:numel(freeP)
    Sp.exp.nor(i)=Sp.exp.raw(1)/Sp.exp.raw(i);
    Sp.es.nor(i)=Sp.es.raw(i)/Sp.es.raw(1);
end
fprintf('\n\n\n\n\n %15s %15s \n','Std.exp.norm','Std.est.norm')
for i=1:numel(Sp.exp.nor)
    fprintf('\n %15i %15f \n',[Sp.exp.nor(i)';Sp.es.nor(i)']);
end

%plot
figure('color',[1 1 1]);
hold all
myerrorbar(Sp.exp.nor,Sp.es.nor,'yError',Sp.std,...
    'Symbol=o',...
    'Markersize=30',...
    'Color=[0 0 0]')
plot(Sp.exp.nor,Sp.exp.nor,'k:');

%graphics
title('Subject representations of prior width')
xlabel('Experimental std of prior (ratio to weaker prior)')
ylabel('Estimated std of prior (ratio to weaker prior)')
ylim([1 10])
xlim([0 max(Sp.exp.nor)+2])
axis square
set(gca,'fontsize',12)
set(gca,'xtick',Sp.exp.nor,'xticklabel',{'80/80','80/40','80/20','80/10'})
drawPublishAxis
function Sp=GetPrior2a(freeP,factor)

%INPUTS
%freeP: In descending order
%factor: factors associated to the freeP

%Convert cell input to matrix
if iscell(freeP)==1
    freeP=[freeP{:}];
end

%Remove cases of infinite prior from the analysis
if freeP(:,1)==+inf
    freeP(1)=[];
end

%Representation of the prior
%.......................................................---
%Get the priors' true and estimated values
Sp.exp.raw     =factor;
Sp.es.raw=freeP(5:8);

%Show summary results
%Variables
fprintf('\n\n\n\n\n %15s %15s \n',...
    'Std.exp.norm',...
    'Std.est.norm')
%Values
for i=1: numel(Sp.exp.raw)
    fprintf('\n %15i %15f \n',...
        [Sp.exp.raw(i)';Sp.es.raw(i)']);
end

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
hold all
title ('Subject representations of prior width')
xlabel('{\sigma}_p_r_i_o_r','fontsize',12)
ylabel('{\sigma}_p_r_i_o_r','fontsize',12)

plot(Sp.exp.raw,Sp.es.raw,'-ko',...
    'markerfacecolor','k',...
    'markersize',15,...
    'displayname','subject' );
plot(Sp.exp.raw,Sp.exp.raw,'k:',...
    'linewidth',6,'displayname','True ratio');

%ylim([1 10])
xlim([0 max(Sp.exp.raw)+2])
lg=legend('location','Northwest');
box(lg,'off')
axis square
set(gca,'fontsize',12)
drawPublishAxis
%%%...........-  TO DO:add a linear fit to subjects' data !!!!!!
%%%...........-  plot also the width of the likelihood !!!!!!!!


%Representation of the likelihood
%.......................................................---
%Get the width of the likelihood
R0   =freeP(1);
Rmax =freeP(2);
n    =freeP(3);
C50  =freeP(4);

%Set factor 1,e.g.,coherence
c=[0.06 0.12 0.24];
%loop over levels of factor 1
for i=1:numel(c)
    sl(i)=1./(R0 + Rmax.*((c(i).^n)./(C50 + (c(i).^n))));
end
figure('color','w');
hold all
plot(c,sl,'ko','markerfacecolor','k','markersize',8,'linewidth',2)
ylabel('Std of the likelihood','fontsize',20)
xlabel('Coherence','fontsize',20)
set(gca,'fontsize',20)
drawPublishAxis
function Sp=GetPrior5(freeP,stdPa,factor,factor2)

%Get the priors' true and estimated values
Sp.exp.raw=factor;
Sp.es.raw=freeP(4:7);
if isempty(stdPa)==0;
    Sp.std=stdPa(4:7);
else
    Sp.std=[];
end

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
title ('Internal prior strength','fontsize',12)
xlabel('Experimental k_p_r_i_o_r','fontsize',12)
ylabel('Estimated k_p_r_i_o_r','fontsize',12)
hold all
myerrorbar(Sp.exp.raw,Sp.es.raw,'yError',Sp.std,...
    'Symbol=o',...
    'Markersize=30',...
    'Color=[0 0 0]');
plot(Sp.exp.raw,Sp.exp.raw,'k:','linewidth',5)
linefit(Sp.exp.raw,Sp.es.raw','k')
%set(gca,'xtick',0:10:90,'xticklabel',0:10:90,'fontsize',20)
drawPublishAxis


%Get the llh' true and estimated values
Sl.es.raw=freeP(1:3);
if isempty(stdPa)==0;
    Sl.std=stdPa(1:3);
else
    Sl.std=zeros(numel(Sl.es.raw),1);
end

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
hold all
title ('Likelihood strength','fontsize',12)
xlabel('Coherence','fontsize',12)
ylabel('Estimated k_l_l_h','fontsize',12)
myerrorbar(factor2,Sl.es.raw,'yError',Sl.std,'Color=[0 0 0]');
xlim([0 .3])
ylim([0 max(Sl.es.raw)+max(Sl.std)])
set(gca,'xtick',0:0.06:.3,'xticklabel',0:0.06:.3,'fontsize',20)
drawPublishAxis


%Support functions
%calculate the circular statistics of the data
function data=statcircular(coord)
%input a vector of cartesian coordinates (coord)

%register the coordinates of the input directions
data.coord.all=coord;

%convert from cartesian coordinates to angles (in degree)
data.deg.all=getangle(coord(:,1),coord(:,2));

%calculate the cartesian coordinates of the mean direction est
data.coord.mean=nanmean(coord,1);

%calculate the mean direction estimate (in degree)
data.deg.mean=getangle(data.coord.mean(:,1),data.coord.mean(:,2));

%calculate the std to the mean direction est (in degree);!!! could be a
%subfunction itself
%Apply the rule of thumb that follows. It seems to work fine intuitively. It would be nice to
%fine a cleaner way to calculate the std.
%initialize the 'sample' and 'mean' variables used to calculate the std
data.num=numel(data.deg.all);%sample size
data.deg.allforstd=data.deg.all;
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

%if the resulting mean direction is between 0 and 180.
if data.deg.mean+180<=360
    %if estimation is>=mean direction+180
    for i=1:data.num
        if data.deg.all(i)>=data.deg.mean+180
            data.deg.allforstd(i)=data.deg.all(i)-360;
        end
    end
    %if the resulting mean direction is between 180 and 360.
else
    %if the estimated direction sampled is <=the mean direction - 180
    for i=1:data.num
        if data.deg.all(i)<=data.deg.mean-180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

%now calculate the variance of the estimated direction.
data.deg.var=nanmean((data.deg.allforstd-data.deg.meanforstd).^2,1);

%and now calculate the std
data.deg.std=sqrt(data.deg.var);

%and now calculate the sem
data.deg.sem=data.deg.std/sqrt(data.num);
%Convert from cartesian coordinates to angles (degree)
function [angle]=getangle(x,y)
% check! to check if the function works fine, write
% e.g., input=180; output=getangle(cos(angle*pi/180),sin(angle*pi/180));
% if the function works input and output should always be the same between
% 0 and 360.

% convert from cartesian coordinates to angle in radians
angle = atan(y./x); % theta = arctan(opposite/adjacent);

% adjust each angle according to his quadrant (in degree)
for i=1:numel(x) % sample each angle
    if x(i)>=0 && y(i)>=0                   %(quadrant 1)
        angle(i) = angle(i)*180/pi;
    elseif x(i)<0                      %(quadrant 2 & 3)
        angle(i) = angle(i)*180/pi + 180;
    elseif x(i)>=0 && y(i)<0               %(quadrant 4)
        angle(i) = angle(i)*180/pi + 360;
    end
end
%Calculate the angle formed by two vectors (distanceR2
function angle=vectors2signedAngle(v1,v2)
%Inputs are two vectors' coordinates
%xV1=v1(1)
%yV1=v1(2)
%xV2=v2(1)
%yV2=v2(2)

%e.g.,v1.x=0;v1.y=1;v2.x=1;v2.y=0;
%angle=- (180/pi) * atan2(v1.x*v2.y - v1.y*v2.x,v1.x*v2.x+v1.y*v2.y)
%gives 90 degrees.

%v1.x=1;v1.y=0;v2.x=0;v2.y=1;
%angle=-(180/pi)*atan2(v1.x*v2.y-v1.y*v2.x,v1.x*v2.x+v1.y*v2.y)
%gives - 90 degrees.

%data
xV1=v1(:,1);
yV1=v1(:,2);
xV2=v2(:,1);
yV2=v2(:,2);

%Calculate the angle in degree separating the two vectors
angle=-(180/pi).*atan2(xV1.*yV2-yV1.*xV2,xV1.*xV2+yV1.*yV2);
%Fit the data with a linear model
function linefit(x,y,linecolor)

%available data
i=~isnan(y);
y=y(i);

%fit
xfit=x;
P=polyfit(x(i),y,1);
yfit=polyval(P,x);
plot(xfit,yfit,'-',...
    'linewidth',0.5,...
    'color',linecolor);
%Plot errorbars
function retval=myerrorbar(x,y,varargin)

%Plot errorbar
%myerrorbar.m
%
%     usage: myerrorbar(x,y,varargin)
%        by: justin gardner
%      date: 06/24/07
%   purpose: draw plots with error bars
%      e.g.: y error bars
%            myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10));
%            x error bars
%            myerrorbar(1:10,rand(1,10),'xError',0.5*rand(1,10));
%            x and yerror bars
%            myerrorbar(1:10,rand(1,10),'yError',0.5*rand(1,10),'xError',0.5*rand(1,10));
%            different lower and upper bounds
%            myerrorbar(1:10,rand(1,10),'yLow',2*rand(1,10),'yHigh',0.5*rand(1,10));
%   options: Symbol=symbol to use,default 'o-'
%            Color=symbol color,default 'k'
%            MarkerFaceColor=symbol face color,defaults to Color
%            MarkerEdgeColor=symbol edge color,defaults to Color
%            MarkerSize=symbol size,default 8
%            myerrorbar(1:10,rand(1,10),'yError',rand(1,10)/2,'Symbol=s-','Color=[1 0.5 0]');
%            tee=draw tees or not,default 0
%            yTeelen=length of tee on y error,default to 1/10 of x spacing
%            xTeelen=length of tee on x error,default to 1/10 of y spacing

%check arguments
if nargin < 2
    help myerrorbar
    return
end

%check for old style usage
if (nargout==1) || ((length(varargin) >=1) && isnumeric(varargin{1}))
    retval=myerrorbarold(x,y,varargin);
    return
end

%get arguments
getArgs(varargin);

%no passed in x
if ieNotDefined('x');x=1:length(y);end

%get length of x
n=length(x);

%get y upper and lower bounds
if ~ieNotDefined('yError'),yLow=yError;yHigh=yError;end
if ieNotDefined('yLow'),yLow=zeros(1,n);end
if ieNotDefined('yHigh'),yHigh=yLow;end
if ieNotDefined('yErrorBarType') yErrorBarType='both';end

%get x upper and lower bounds
if ~ieNotDefined('xError'),xLow=xError;xHigh=xError;end
if ieNotDefined('xLow'),xLow=zeros(1,n);end
if ieNotDefined('xHigh'),xHigh=xLow;end

%colors and symbols
if ieNotDefined('Symbol'),Symbol='o-';end
if ieNotDefined('Color')
    if ~ieNotDefined('MarkerFaceColor')
        Color=MarkerFaceColor;
    else
        Color='k';
    end
end
if ieNotDefined('MarkerEdgeColor'),MarkerEdgeColor=Color;end
if ieNotDefined('MarkerFaceColor'),MarkerFaceColor=Color;end
if ieNotDefined('MarkerSize'),MarkerSize=8;end
if ieNotDefined('LineWidth'),LineWidth=0.5;end

%whether to draw tees or not
if ieNotDefined('tee'),tee=0;end
if tee
    if ieNotDefined('yTeelen'),yTeelen=mean(diff(x))/10;end
    if ieNotDefined('xTeelen'),xTeelen=mean(diff(y))/10;end
end
hold on
%plot the y error bars
if any(any(yLow ~=0)) || any(any(yHigh ~=0))
    for i=1:length(x)
        switch yErrorBarType
            case {'both','b'}
                plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
            case {'lower','lo','l','bottom','bot'}
                plot([x(i) x(i)],[y(i)-yLow(i) y(i)],'-','Color',Color,'LineWidth',LineWidth);
            case {'higher','upper','up','top','hi'}
                plot([x(i) x(i)],[y(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
            case {'logy'}
                if (y(i)-yLow(i)) > 0
                    plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
                else
                    disp(sprintf('(myerrorbar) Dropping lower errorbar on %i which goes to %f',i,y(i)-yLow(i)));
                    plot([x(i) x(i)],[y(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
                end
            otherwise
                plot([x(i) x(i)],[y(i)-yLow(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
        end
        %draw the tees if necessary
        if tee
            plot([x(i)-yTeelen/2 x(i)+yTeelen/2],[y(i)-yLow(i) y(i)-yLow(i)],'-','Color',Color,'LineWidth',LineWidth);
            plot([x(i)-yTeelen/2 x(i)+yTeelen/2],[y(i)+yHigh(i) y(i)+yHigh(i)],'-','Color',Color,'LineWidth',LineWidth);
        end
    end
end

%plot the x error bars
if any(any(xLow ~=0)) || any(any(xHigh ~=0))
    for i=1:length(x)
        plot([x(i)-xLow(i) x(i)+xHigh(i)],[y(i) y(i)],'-','Color',Color,'LineWidth',LineWidth);
        %draw the tees if necessary
        if tee
            plot([x(i)-xLow(i) x(i)-xLow(i)],[y(i)-xTeelen/2 y(i)+xTeelen/2],'-','Color',Color,'LineWidth',LineWidth);
            plot([x(i)+xHigh(i) x(i)+xHigh(i)],[y(i)-xTeelen/2 y(i)+xTeelen/2],'-','Color',Color,'LineWidth',LineWidth);
        end
    end
end

%plot the symbols
plot(x,y,Symbol,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'Color',Color,'MarkerSize',MarkerSize,'LineWidth',LineWidth);
%axis
function axs=makeAxes(F)
%axes width
width=1/(F.g1.numL+1);

%space between axes
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);

%positionne axes
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end
%von Mises pdf
function mPdfs=vmPdfs(x,u,k)
%k must be in the range [0:709]. when k>709, probabilities are inf.
%k must be in the range [0:700]. when k>700, besseli is inf amd mPdfs is 0.
%When normalizing the von mises to probabilities we get NaN values
%(0/sum(0)).
w=1;

%radians
x=de2r(x,1); x=x';
u=de2r(u,1); u=u';
x2=x(:,ones(numel(u),1));
u2=u(ones(numel(x),1),:);
k2=k(ones(numel(x),1),:);
mPdfs=exp(k2.*cos(w*(x2-u2)))./(2*pi.*besseli(0,k2));
%pdf
Z_=sum(mPdfs);
Z=Z_(ones(numel(x),1),:);
mPdfs=mPdfs./Z;
%vector average
function [deg,rad,cart]=vectorAverage(v)
%v is a directional vector (cartesian coordinates)
%r is the unit circle radius.

% %cartesian of the mean
% Av.cart=nanmean(v,1);
% x=Av.cart(:,1);
% y=Av.cart(:,2);
%
% %radian (arctan(opposite/adjacent))
% %Av.rad=atan(y./x);
% [Av.rad, Av.deg]=cart2polar(Av.cart);

%cartesian of the mean
cart=nanmean(v,1);
x=cart(:,1);
y=cart(:,2);

%radian (arctan(opposite/adjacent))
%Av.rad=atan(y./x);
[rad, deg]=cart2polar(cart);
%polar to cartesian
function v=polar2cartesian(di,r)
%di:angs in degree
%r:radius of the unit circle

%degrees & radians
di2.deg=di;
di2.rad=(di2.deg*pi)/180;
% di2.rad=(di2.deg*3.14)/180;


%cartesian vectors
x=r.*cos(di2.rad);
y=r.*sin(di2.rad);
v=[x y];
%polar to cartesian
function [rad,deg]=cart2polar(v)
%v is the cartesian vector

%radian (arctan(opposite/adjacent))
x=v(:,1);
y=v(:,2);
rad=atan(y./x);

%degrees
deg=ra2d(rad);

%adjust each angle according to his quadrant (in degree) otherwise atan
%gives wrong values for negative x and y.
for i=1:numel(x)
    if x(i)>=0 && y(i)>=0              %(quadrant 1)
        deg(i)=deg(i);
    elseif x(i)<0                      %(quadrant 2 & 3)
        deg(i)=deg(i)+180;
    elseif x(i)>=0 && y(i)<0           %(quadrant 4)
        deg(i)=deg(i)+360;
    end
end
%degrees to radians
function radians=de2r(ang,sign)

%not signed radians (1:2*pi)
radians=(ang/360)*2*pi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
end
%convert radians to degrees
function degrees=ra2d(theta)

%When input radians are between 0:2*pi
degrees=(theta/(2*pi))*360;

%if larger than 360 degrees then subtract
%360 degrees
while (sum(degrees>360))
    degrees = degrees - (degrees>360)*360;
end

%if less than 360 degreees then add
%360 degrees
while (sum(degrees<-360))
    degrees = degrees + (degrees<-360)*360;
end

%When radians are signed between -pi:pi.
degrees(degrees<0)=degrees(degrees<0)+360;
%Girshick lookup table
function ul=lookuptable(ulall,up,kl,kp)
%method2:row upo/col ul; fun is the function to minimize to find ul
%solution(degrees).ulall2, up and upo must be in signed radians because the von
%mises from which the equation is derived takes radians values for ulall2,
%and up.
%lookuptable inputs are in degrees and are converted in radians
% upo=[1:0.1:360]';
upo=[1:1:360]';

upo2=upo(:,ones(1,numel(ulall)));
ulall2=ulall(ones(numel(upo),1),:);
fun=abs(de2r(upo2,1)-(de2r(ulall2,1)+atan2(sin(de2r(up,1)-de2r(ulall2,1)),(kl/kp)+cos(de2r(up,1)-de2r(ulall2,1)))));
[dummy,I]=min(fun,[],2);

%in degrees
ul=ulall(I);