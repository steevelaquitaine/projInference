

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



function [Sp,pred,fitP,R2,udata,sdata,dataOut,predOut,Fbp,F,stdPa,fitPt]=modelBayesInf4(data,Fs,fig,FofLl2del,Ll2del)

%.......
%OUTPUTS (clean up!)
%.......
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

%Organize databank
fprintf('\n %12s \n','Now create data matrix...')
Makedatabank(data,FofLl2del,Ll2del);

%Get factors
fprintf('\n %12s \n','Now sorting the factors...')
[datafitm,datafitt,F]=initFs(Fs,fig);

%........
%Modeling
%........

%model 1

%Fit data's mean with width of llh(std) an inverse function of coherence.
%You are free to set (or not,[]) 1 free parameters.
%fitP=[0.1 0.2 0.3 0.4];
%[pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut]=LSft(datafitm,fitP);%Fit
%drawPred(pred,udata,Fbp,F,fig);%Plot model predictions
%Sp=GetPrior1(fitP,F.g2.L);%Get prior's representation
%----- Try to fit also to std of the data here -----!


%model 2a

%Fit data's mean and std with std of llh(std) a hyperbolic function of
%coherence.
%    4 free parameters for hyperbolic function,std of priors and motor noise
%display('Now fitting the model....')
%You are free to set (or not,[]) 9 free parameters.
%fitP=[0.1 0.2 0.3 0.4];
%[pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut]=LSft2a(datafitm,fitP);
%display('Now drawing the predictions....')
%drawPred5(pred,udata,sdata,Fbp,F,fig)
%Sp=GetPrior2a(fitP,F.g2.L);%Get prior's representation


%model 3

%Fit data's mean with width of llh(std) as a free parameter and assuming the true priors
%[pred,udata,Fbp,fitP,R2]=LSft3(datafitm);
%----- Try to fit also to std of the data here -----!



%model 4

%Fit data's mean and std with width of llh(std) as a free parameter and assuming the true priors
%[pred,udata,sdata,Fbp,fitP,R2]=LSft4(datafitm);
%[fig1]=drawPred4(pred,udata,sdata,Fbp,F,fig);


%model 5 (~30 min)

%Fit data's mean and std with 
    %std of llh,std of priors and motor noise as 8 free parameters
%fprintf(' \n %s \n','Now fitting the model with the mean and std of the data ...')
%You are free to set (or not,[]) your own free parameters.
%fitP=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.01];
%[pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa]=LSft5(datafitm,fitP);
%fprintf(' \n %s \n','Now drawing model predictions....')
%fig1=drawPred5(pred,udata,sdata,Fbp,F,fig);
%Sp=GetPrior5(fitP,stdPa,F.g2.L,F.g1.L);


%model 6 (~48 min)

%Fit mean and std with Bayesian inference
    %std llh,std priors and as 7 free parameters
    % there is no motor noise because I don't know yet how to implement it.
fprintf('\n %s \n','Now fitting the model with the mean and std of the data ...')
fitP=[];%[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.01];
[pred,data,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa,fitPt]=LSft6(datafitt,fitP,F);
fprintf('\n %s \n','Now drawing model predictions....')
fig1=drawPred6(pred,F,fig,datafitt);

%...............
%Create databank
%...............
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

%.............
%Sort factors)
%.............
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

%............
%Process data
%............
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

%........
%Fit data
%........
function [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut]=LSft(datafitm,fitP)
%GOAL
    %Fit mean with one free parameter
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
datafitm.data(cellfun(@(l) any(isnan(l)),l),:)=[];

%Collect data
udata=[datafitm.data{:,1}]';
sdata=[datafitm.data{:,2}]';

%Factors
%g1
F.g1=[datafitm.data{:,3}]';
F.g1L=unique(F.g1);
F.g1L=sort(F.g1L,'descend');%order

%g2
F.g2=[datafitm.data{:,4}]';
F.g2L=unique(F.g2);
F.g2L=sort(F.g2L,'descend');%order

%g3
F.g3=[datafitm.data{:,5}]';
F.g3L=unique(F.g3);
F.g3L=sort(F.g3L,'ascend');%order

%Store
Fbp=[F.g1 F.g2 F.g3];


%If there is no parameters,fit the model to the data
if isempty(fitP)==1
    
    %Initial parameters
    k0_p1=0.05:0.05:0.15;
    k0_p2=0.05:0.05:0.15;
    k0_p3=0.05:0.05:0.15;
    k0_p4=0.05:0.05:0.15;
    
    %Fit
    %.........................................................................................................................--
    %iteration=0;
    %Loop over free parameters
    for i=1:numel(k0_p1)%e.g.,prior 1
        for j=1:numel(k0_p2)%e.g.,prior 2
            for k=1:numel(k0_p3)%e.g.,prior 3
                for l=1:numel(k0_p4)%e.g.,prior 4
                    
                    %Fit
                    fitPtmp=fmincon( @(fitPtmp) makeSSE1(udata,Fbp,...
                        fitPtmp),...
                        [k0_p1(i);k0_p2(j);k0_p3(k);k0_p4(l)],...
                        [],[],[],[],...
                        [0 0 0 0],...
                        [+inf +inf +inf +inf],...
                        []);
                    
                    %Store
                    fitPbkp{i,j,k,l}=fitPtmp;
                    
                    %Calculate SSE
                    SSE_bkp(i,j,k,l)=makeSSE1(udata,Fbp,fitPtmp);
                    
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
    R2=makeR22(udata,pred00.mean');
    
    %Store output
    dataOut=udata;
    predOut=pred00.mean';
    
    %Get model predictions on a larger space
    pred=makePre1(F,fitP);
    
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
datafitm.data(cellfun(@(l) any(isnan(l)),l),:)=[];

%Collect data
udata=[datafitm.data{:,1}]';
sdata=[datafitm.data{:,2}]';

%Factors
%g1
F.g1=[datafitm.data{:,3}]';
F.g1L=unique(F.g1);
F.g1L=sort(F.g1L,'descend');%order

%g2
F.g2=[datafitm.data{:,4}]';
F.g2L=unique(F.g2);
F.g2L=sort(F.g2L,'descend');%order

%g3
F.g3=[datafitm.data{:,5}]';
F.g3L=unique(F.g3);
F.g3L=sort(F.g3L,'ascend');%order

%Store
Fbp=[F.g1 F.g2 F.g3];

%If there is no parameters,fit the model to the data
if isempty(fitP)==1
    
    %Initial parameters
    R0_0=0%:0.05:0.1;
    Rmax_0=1%:0.05:0.1;
    n_0=1%:0.05:0.1;
    C50_0=1%:0.05:0.1;
    sp0_1=1:81:164;
    sp0_2=1:81:164;
    sp0_3=1:81:164;
    sp0_4=1:81:164;
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
                                        [0 0 0 0 0 0 0 0 0],...
                                        [+inf +inf +inf +inf +inf +inf +inf +inf +inf],...
                                        []);
                                    
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
    R2=makeR22(udata,pred00.mean');
    
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
datafitm.data(cellfun(@(l) any(isnan(l)),...
    l),:)=[];


%Fit sl for each coherence
%for data(coh==1),set the priors as true and fit sl.

%Collect the data for coherence ii
Ftofit.labels=cell2mat(datafitm.data(:,3));%(e.g.,coh)
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
    udata{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii),1));
    
    %Conditions
    %g1
    F.g1{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii) ,3));
    %g2
    F.g2{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii) ,4));
    %g3
    F.g3{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii) ,5));
    %Store
    Fbp{ii}=[F.g1{ii}  F.g2{ii}  F.g3{ii}];
   
    %.........................................................................................................................--
    %Fit data
    %.........................................................................................................................--
    %Loop over the initial parameters
    for i=1:numel(k0)
        
        %Fitting
        fitPtmp=fmincon( @(fitPtmp) makeSSE3(udata{ii},Fbp{ii},...
            fitPtmp),...
            k0(i),...
            [],[],[],[],0,+inf,[]);
        
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
datafitm.data(cellfun(@(l) any(isnan(l)),...
    l),:)=[];

%Fit sl for each coherence
%for data(coh==1),set the priors as true and fit sl.

%Collect the data for coherence ii
Ftofit.labels=cell2mat(datafitm.data(:,3));%(e.g.,coh)
Ftofit.L=sort(unique(Ftofit.labels),'descend');

%Initial parameters
k0=0.05%:0.05:20;


fitPBkp=[];
%Loop over coherences
for ii=1:numel(Ftofit.L);

    %.........................................................................................................................--
    %Select variables
    %.........................................................................................................................--
    %Data
    udata{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii),1));
    sdata{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii),2));

    %Conditions
    %g1
    F.g1{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii) ,3));%e.g.,coh
    %g2
    F.g2{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii) ,4));
    %g3
    F.g3{ii}=cell2mat(datafitm.data(Ftofit.labels==Ftofit.L(ii) ,5));
    %Store
    Fbp{ii}=[F.g1{ii}  F.g2{ii}  F.g3{ii}];
   
    %.........................................................................................................................--
    %Fit data
    %.........................................................................................................................--
    %Loop over the initial parameters
    for i=1:numel(k0)
        
        %Fitting
        fitPtmp=fmincon( @(fitPtmp) makeSSE4(udata{ii},sdata{ii},Fbp{ii},...
            fitPtmp),...
            k0(i),...
            [],[],[],[],0,+inf,[]);
        
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
function [pred,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa]=LSft5(datafitm,fitP)
%Fit mean and variance with std(llh) as a free parameter and assuming the true priors
%Parameters search (model fitting by Least Square Optimization (LSO))

%Remove "NaN" data
l=datafitm(:,1);
datafitm.data(cellfun(@(l) any(isnan(l)),l),:)=[];

%.........................................................................................................................--
%Select variables
%.........................................................................................................................--
%Collect data
udata=[datafitm.data{:,1}]';
sdata=[datafitm.data{:,2}]';

%Conditions
%g1
F.g1=[datafitm.data{:,3}]';
F.g1L=unique(F.g1);
F.g1L=sort(F.g1L,'descend');%order

%g2
F.g2=[datafitm.data{:,4}]';
F.g2L=unique(F.g2);
F.g2L=sort(F.g2L,'descend');%order

%g3
F.g3=[datafitm.data{:,5}]';
F.g3L=unique(F.g3);
F.g3L=sort(F.g3L,'ascend');%order

%Store
Fbp=[F.g1 F.g2 F.g3];


%If there is no parameters,fit the model to the data
if isempty(fitP)==1
    
    %Initial parameters(Elapsed time is 2114.179696 seconds)
    sl1_0=1:81:164;%std of likelihood: c=0.24
    sl2_0=1:81:164;
    sl3_0=1:81:164;
    sp1_0=1:81:164;%std of prior: 80 deg
    sp2_0=1:81:164;
    sp3_0=1:81:164;
    sp4_0=1:81:164;
    sM_0=0;%1:81:164;%motor noise
    
    %.........................................................................................................................--
    %Fitting
    %.........................................................................................................................--
    %iteration=0;
    %Loop over free parameters
    for i=1:numel(sl1_0)
        for j=1:numel(sl2_0)
            for k=1:numel(sl3_0)
                for l=1:numel(sp1_0)
                    for m=1:numel(sp2_0)
                        for n=1:numel(sp3_0)
                            for e=1:numel(sp4_0)
                                for f=1:numel(sM_0)
                                    
                                    %                                    %Fit
                                    %                                    [fitPtmp]=fmincon( @(fitPtmp) makeSSE5(udata,sdata,Fbp,...
                                    %                                        fitPtmp),...
                                    %                                        [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e);sM_0(f)],...
                                    %                                        [],[],[],[],...
                                    %                                        [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001],...
                                    %                                        [+inf +inf +inf +inf +inf +inf +inf +inf],[]);
                                    %
                                    %Fit
                                    [fitPtmp,~,~,~,~,~,Hessian]=fmincon( @(fitPtmp) makeSSE5(udata,sdata,Fbp,...
                                        fitPtmp),...
                                        [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e);sM_0(f)],...
                                        [],[],[],[],...
                                        [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001],...
                                        [+inf +inf +inf +inf +inf +inf +inf +inf],[]);
                                                                        
                                    %Store
                                    fitPbkp{i,j,k,l,m,n,e,f}=fitPtmp;
                                    
                                    %Calculate the SSE
                                    SSE_bkp(i,j,k,l,m,n,e,f)=makeSSE5(udata,sdata,Fbp,fitPtmp);
                                    
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
    
    
    %Get the lower SSE
    [minSSE,position]=min(SSE_bkp(:));
    [i,j,k,l,m,n,e,f]=ind2sub(size(SSE_bkp),position);
    
    %Get the best parameters
    fitP=fitPbkp{i,j,k,l,m,n,e,f};
    
    %Get extended model's predictions
    [pred,~]=makePre(F,fitP);
    
    %Calculate data-restricted model predictions
    [~,pred00]=makeSSE5(udata,sdata,Fbp,fitP);
    
    %Get the R^2
    %%method 1
    R2=makeR2([udata;sdata],minSSE);
    %%method 2
    %[~,pred2]=makeSSE5(udata,sdata,Fbp,fitP);
    %R2=makeR22([udata;sdata],[pred2.mean';pred2.std']);
    
    %Store output
    dataOut=[udata;sdata];
    predOut=[pred00.mean';pred00.std'];
    
    %Get standard deviation of model parameters
    stdPa=GetStdP(Hessian,dataOut,predOut,fitP);
    
%If free parameters are input,only draw predictions
elseif isempty(fitP)==0
    
    %Calculate the SSE
    [~,pred00]=makeSSE5(udata,sdata,Fbp,fitP);
    
    %Calculate the R2
    R2=makeR22([udata;sdata],[pred00.mean';pred00.std']);
    
    %Store output
    dataOut=[udata;sdata];
    predOut=[pred00.mean';pred00.std'];
    
    %Get the model's predictions generalized
    [pred,~]=makePre(F,fitP);
    
    %Getting standard deviation of model parameters is not meaningful here
    stdPa=[];
end
function [pred,data,udata,sdata,Fbp,fitP,R2,dataOut,predOut,stdPa,fitPt]=LSft6(datafitt,fitP,F)
%This function fit trial data from "datafitt" with bayesian inference

%get data and factors
%check normalization
if ~isfield(F.g3,'distthisT')
    data=datafitt.data{strcmp(datafitt.nm,'dir')==1};
    Fbp=[F.g1.thisT F.g2.thisT F.g3.thisT];
    Pmean=225;
else
    data=datafitt.data{strcmp(datafitt.nm,'distances')==1};
    Fbp=[F.g1.thisT F.g2.thisT F.g3.distthisT'];
    Pmean=0;
end
%mean and std
udata=datafitt.data{strcmp(datafitt.nm,'mean')==1};
sdata=datafitt.data{strcmp(datafitt.nm,'std')==1};


%when no parameters are input,fit.
if isempty(fitP)==1
    
    %Options: fmincon complains with the default soLer 
    %('Trust-reflective-region') and switch to active-set. 
    %So I directly use active-set.
    options=optimset('Display','off','MaxFunEvals',50000,...
                     'TolFun',1e-06,...
                     'MaxIter',50000,...
                     'Algorithm','active-set');
    
    %Initial fit parameters
    %llhs std
    sl1_0=1:81:164;
    sl2_0=1:81:164;
    sl3_0=1:81:164;
    %priors std
    sp1_0=1:81:164;
    sp2_0=1:81:164;
    sp3_0=1:81:164;
    sp4_0=1:81:164;
    %motor noise
    sM_0=0;
    
    %Fitting
    %iteration=0;
    %Loop over fit parameters.
    for i=1:numel(sl1_0)
        for j=1:numel(sl2_0)
            for k=1:numel(sl3_0)
                for l=1:numel(sp1_0)
                    for m=1:numel(sp2_0)
                        for n=1:numel(sp3_0)
                            for e=1:numel(sp4_0)
                                for f=1:numel(sM_0)
                                    
                                    %Fit mean and std.
                                    [fitPtmp,~,~,~,~,~,Hessian]=fmincon(@(fitPtmp) makeSSE6(data,Fbp,fitPtmp,Pmean),...%fun
                                        [sl1_0(i);sl2_0(j);sl3_0(k);sp1_0(l);sp2_0(m);sp3_0(n);sp4_0(e);sM_0(f)],...%x0
                                        [],[],[],[],...%Aeq,Beq
                                        [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001],...%lb
                                        [+10000 +10000 +10000 +10000 +10000 +10000 +10000 +10000],...%ub
                                        [],...%nonlcon
                                        options);
                                                                        
                                    %Store fit parameters
                                    fitPbkp{i,j,k,l,m,n,e,f}=fitPtmp;
                                    
                                    %Calculate the SSE
                                    SSE_bkp(i,j,k,l,m,n,e,f)=makeSSE6(data,Fbp,fitPtmp,Pmean);
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
        toc
    end
    
    %Get the lower SSE
    [minSSE,position]=min(SSE_bkp(:));
    [i,j,k,l,m,n,e,f]=ind2sub(size(SSE_bkp),position);
    
    %Get the best parameters
    fitP=fitPbkp{i,j,k,l,m,n,e,f};
    
    %Calculate data-restricted model predictions
    [~,pred00,~,fitPt]=makeSSE6(data,Fbp,fitP,Pmean);
    
    %Get general model's predictions
    [pred,~]=makePre6(F,fitPt,Pmean);
    
    %Get R^2
    %method 1
    R2=makeR2(data,minSSE);
    %method 2
    %[~,pred2]=makeSSE5(data,sdata,Fbp,fitP);
    %R2=makeR22([data;sdata],[pred2.mean';pred2.std']);
    
    %Store output
    dataOut=data;
    predOut=pred00.mean;
    
    %Get standard deviation of model parameters
    stdPa=GetStdP(Hessian,dataOut,predOut,fitP);
    
%If fit parameters are input,just draw predictions
elseif isempty(fitP)==0
    
    %Calculate the SSE
    [~,pred00,~,fitPt]=makeSSE6(data,Fbp,fitP,Pmean);
    
    %Calculate the R2
    R2=makeR22(data,pred00.mean);
    
    %Store output
    dataOut=data;
    predOut=pred00.mean;
    
    %Get the model's predictions generalized
    [pred,~]=makePre6(F,fitPt,Pmean);
    
    %std of fit parameters have no meaning here.
    stdPa=[];
end

%................
%Draw predictions
%................
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
    p10=plot(F.g3.L,F.g3.L','k:',...
        'linewidth',2,...
        'Displayname','Ideal predictions');
    
    %organize data for plotting
    %.........................................................................................................................---
    %change the level of factor 2 (e.g.,prior's strength)
    for j=1:F.g2.numL
        %change the level of factor 3 (displayed directions)
        %coordinates of single conditions
        indX.g1g2(j,k)={intersect( indX.g1.Ll_i{k},indX.g2.Ll_i{j} ) };
        %store data
        fig1.datamean{j,k}=udata(indX.g1g2{j,k});%es{i,j,k}.deg.mean;

        %extract information about data
        %g1
        fig1.dataInfog1{j,k}=Fbp(indX.g1g2{j,k},1);%F.g1.L(k);
        %g2
        fig1.dataInfog2{j,k}=Fbp(indX.g1g2{j,k},2);%F.g2.L(j);
        %g3
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
        p13a=plot(F.g3.L,pred.mean(:,j,k),'-',... %groups 3 forms x-axis
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',2,...
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
% F.g2.color={[0.5 0 0],...
%     c(55,:),...
%     [1 0.4 0],...
%     [0.75 0.75 0],...
%     c(40,:),...
%     c(32,:),...
%     c(27,:),...
%     c(22,:),...
%     c(8 ,:),...
%     c(5 ,:),...
%     c(1 ,:)};%group 2(e.g.,coherences)

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
function [fig1]=drawPred5(pred,udata,sdata,Fbp,F,fig)
global databank

%Draw group1-Ll1 (e.g.,prior std=80)
%Enlarge figure for good quality publication
%fig1.hdle=figure('Position',[0 0 1000 400]);%pixels
fig1.hdle=figure('color',[1 1 1]);
fig1.nm=[fig.nm,'_DataAndModel'];


%...........!!!!!! Seriously check the correspondence of udata,sdata and
%Fbp !!!!!!!!!!!!!!!!!!!!!!!!



%.........................................................................................................................
%Set the factors
%.........................................................................................................................
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


%.........................................................................................................................----
%Plot 
%.........................................................................................................................----
%Set the graph's parameters
%colors
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

%if not enough colors
if numel(F.g2.color) < F.g2.numL
    disp(['--- You may want to add ',...
        num2str(F.g2.numL - numel(F.g2.color)),...
        'more colors to the color code ---']);
    return
end

%axes' positions
width=1/(F.g1.numL+1);
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end

%space axes
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end


%Mean

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
    
    %Draw the predictions of ideal performances
    p10=plot(F.g3.L,F.g3.L','k:',...
        'linewidth',1,...
        'Displayname','Ideal predictions');
    
    %organize data for plotting
    %...................................................
    %loop over factor 2's levels (e.g.,prior's strength)
    for j=1:F.g2.numL
        
        %Loop over factor 3's levels (displayed directions)
        %get the positions of each condition
        indX.g1g2(j,k)={intersect( indX.g1.Ll_i{k},indX.g2.Ll_i{j} ) };
        %Store data
        fig1.datamean{j,k}=udata(indX.g1g2{j,k});%es{i,j,k}.deg.mean;

        %extract information about data
        %g1
        fig1.dataInfog1{j,k}=Fbp(indX.g1g2{j,k},1);%F.g1.L(k);
        %g2
        fig1.dataInfog2{j,k}=Fbp(indX.g1g2{j,k},2);%F.g2.L(j);
        %g3
        fig1.dataInfog3{j,k}=Fbp(indX.g1g2{j,k},3);%F.g3.L(i);
              
%        %Normalize data for plotting (because circular data) -already
%        done earlier in the code. should be suppressed
%        %.............................................................
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
        %..........
        p12_=plot(fig1.dataInfog3{j,k},fig1.datamean{j,k} ,'o',... %groups 3 forms x-axis
            'color',F.g2.color{j},...
            'MarkerEdgeColor','w',...
            'markerfacecolor',F.g2.color{j},...
            'markersize',8,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));       
        
        %plot model's predictions
        %.........................
        %mean
        p13a=plot(F.g3.L,pred.mean(:,j,k),'-',... %groups 3 forms x-axis
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',2,...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
        
        %Store data and predictions plotted here....
        
    end

    %set the graph's parameters
    %set the unit steps of the x axis
    xunit=1: 6: F.g3.numL;
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
    xpos(k)=axs.position(k,1) + axs.position(k,3)/2;
    ypos(k)=0.05 + axs.position(k,2)+axs.position(k,4);
    text(xpos(k),...
        ypos(k),...
        strcat(F.g1.nm,': ',num2str(F.g1.L(k))),...
        'fontweight','Bold',...
        'fontsize',12);
end

  



%Plot data - std
%...............

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
    %................
    hold all
    %draw priors' mean
    p1=plot([databank.data{1,7} databank.data{1,7}],[0 150],...
        'b:',...
        'linewidth',1.00005,...
        'DisplayName','Prior mean');
    
    %organize data for plotting
    %..........................
    %Loop over factor 2's levels
    for j=1:F.g2.numL
        
        %Loop over factor 3's levels
        %Get the position of each conditions
        indX.g1g2(j,k)={intersect( indX.g1.Ll_i{k},indX.g2.Ll_i{j} ) };
        
        %Collect the data
        fig1.datastd{j,k}=sdata(indX.g1g2{j,k});
        
        %Extract the factors/levels label
        %g1
        fig1.dataInfog1{j,k}=Fbp(indX.g1g2{j,k},1);%F.g1.L(k);
        %g2
        fig1.dataInfog2{j,k}=Fbp(indX.g1g2{j,k},2);%F.g2.L(j);
        %g3
        fig1.dataInfog3{j,k}=Fbp(indX.g1g2{j,k},3);%F.g3.L(i);
        
%        %Normalize data for plotting (because circular data) -already
%        done earlier in the code. should be suppressed
%        %............................................................
%        %In case the direction is displayed in the 3rd quarter.
%        if fig1.dataInfog3{j,k} > 270 & fig1.datamean{j,k} < 90
%            %"Linearize" the value of the estimated direction
%            fig1.datamean{j,k}=360 + fig1.datamean{j,k};
%            %In case the direction is displayed in the 1st quarter.
%        elseif fig1.dataInfog3{j,k} < 90 & fig1.datamean{j,k} > 270
%            %"Linearize" the value of the estimated direction
%            fig1.datamean{j,k}=fig1.datamean{j,k} - 360;
%        end

       
        %plot the data
        %.............
        p12_=plot(fig1.dataInfog3{j,k},fig1.datastd{j,k},'o',... %groups 3 forms x-axis
            'color',F.g2.color{j},...
            'markerfacecolor',F.g2.color{j},...
            'markersize',8,...
            'MarkerEdgeColor','w',...
            'displayname',strcat(F.g2.nm,':',num2str(F.g2.L(j))));
                
        %plot the model's predictions
        %............................
        %mean
        p13a=plot(F.g3.L,pred.std(:,j,k)','-',... %groups 3 forms x-axis
            'color',F.g2.color{j} - [0.2 0 0],...
            'linewidth',2,...
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
        plot(F.g3.disp,pred.mean(:,j,k),'-',... %groups 3 forms x-axis
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
        p13a=plot(F.g3.disp,pred.std(:,j,k)','-',...
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

%.............
%calculate SSE
%.............
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
SSE=sum((udata-uPo').^2);
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
SSE=sum( ([udata;sdata] - [uPo';sEs']).^2 );

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
SSE=sum( (udata - uPo).^2 );

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
SSE=sum(([udata;sdata]-[uPo;sPo]).^2);
pred.mean=uPo;
pred.std=sPo;
fitP=sl;
function [SSE,pred,fitP]=makeSSE5(udata,sdata,Fbp,fitP)

%...............
%Representation (what)
%...............
%Llh
%mean
uLl=Fbp(:,3);
%std
sl=nan(numel(udata),1);
sl(Fbp(:,1)==0.24)=fitP(1);
sl(Fbp(:,1)==0.12)=fitP(2);
sl(Fbp(:,1)==0.06)=fitP(3);

%Prior
%mean
uPr=225;
%std
sp=nan(numel(udata),1);
sp(Fbp(:,2)==80)=fitP(4);
sp(Fbp(:,2)==40)=fitP(5);
sp(Fbp(:,2)==20)=fitP(6);
sp(Fbp(:,2)==10)=fitP(7);

%Motor noise
sM=fitP(8);

%..........
%Inference (How)
%..........
%mean
uPo=(1./(1+(sl./sp).^2)).*uLl+(1./(1+(sp./sl).^2)).*uPr;
%std
sPo=sqrt(1./((1./sl.^2)+(1./sp.^2)));

%...........
%production
%...........
sEs=sPo + sM;

%...............................
%SSE,predictions and parameters
%...............................
SSE=sum(([udata;sdata]-[uPo;sEs]).^2);
pred.mean=uPo;
pred.std=sEs;%model 6
function [SSE,pred,fitP,fitPt]=makeSSE6(data,Fbp,fitP,fixP)

%...............
%Representation (what)
%...............
%Llh
%mean and std
uLl=Fbp(:,3);
%"see"

sl=nan(numel(data),1);
sl(Fbp(:,1)==0.24)=fitP(1);
sl(Fbp(:,1)==0.12)=fitP(2);
sl(Fbp(:,1)==0.06)=fitP(3);
%"no see"
slnosee=10000;

%Prior
%mean and std
uPr=fixP;
sp=nan(numel(data),1);
sp(Fbp(:,2)==80)=fitP(4);
sp(Fbp(:,2)==40)=fitP(5);
sp(Fbp(:,2)==20)=fitP(6);
sp(Fbp(:,2)==10)=fitP(7);

%Motor noise
sM=fitP(8);

%..........
%Inference (How)
%..........
%first round to get the best llh at each trial.
%Infer with a gaussian llh which std is fixed over trials and with a
%uniform llh at each trial. Compare the two cases' SSE at each trial and
%keep the lower SSE. The description of llh that yield the lower 
%SSE at a given trial explains best subjects'estimation in that trial.
%"see" case
uPosee=(1./(1+(sl./sp).^2)).*uLl+(1./(1+(sp./sl).^2)).*uPr;
%"no see" case
uPonosee=(1./(1+(slnosee./sp).^2)).*uLl+(1./(1+(sp./slnosee).^2)).*uPr;
%get lower SSE llh
sl((data-uPonosee).^2<(data-uPosee).^2)=slnosee;

%make predictions
%percept
uPo=(1./(1+(sl./sp).^2)).*uLl+(1./(1+(sp./sl).^2)).*uPr;
%variability
%sPo=sqrt(1./((1./sl.^2)+(1./sp.^2)));
%production
%sEs=sPo+sM;

%...............................
%SSE,predictions and parameters
%...............................
SSE=sum((data-uPo).^2);
pred.mean=uPo;
%pred.std=sEs;
fitPt=[sl sp repmat(sM,numel(data),1)];

%................
%make predictions
%................
function [pred,F]=makePre1(F,BestfitP)
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

%Model fixed parameters
c =F.g1L;%coherence
uPr=225;%mean of the prior
uLl=F.g3L;

%Model free parameters
k=BestfitP;

%.........................................................................................................................
%Operations
%.........................................................................................................................-
%# in theory,s/c is the variance of the likelihood (i.e.,noise) and can't
%be negative
%# in theory sPr is the variance of the prior and can't be negative.
%# What about k=s/sPr ?

%sPr cannot < 0,thus if k<0,it means s<0. If s<0,s/c>0 only if c<0;
%Thus,when k<0,encoding of coherence<0;
%Not sure how that makes sense....

%Assign the free parameters to each condition
%loop over 
for i=1:numel(F.g1L) %e.g.,coherence 
    for j=1:numel(F.g2L) %e.g.,prior
        for l=1:numel(F.g3L) %e.g.,displayed directions
            
            %Run Bayesian inference (1st hypothesis,mean of the data)
            uPo(l,i,j)=( 1./(1+(k(j)./c(i)).^2) ).*uLl(l)  +  uPr.*( 1./(1+(c(i)./k(j)).^2) );
            
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
function [pred,F]=makePre6(F,BestfitP,fixP)

%single-trial prior
%mean
uPr=fixP;
sp=BestfitP(:,2);

%single-trial llh
%mean
if ~isfield(F.g3,'distL')
    uLl=F.g3.thisT;
else
    uLl=F.g3.distthisT';
end
%std
sl=BestfitP(:,1);

%mean estimate's motor noise
sM=BestfitP(1,3);

%Bayesian infer single-trial estimates
uPo.thisT=(1./(1+(sl./sp).^2)).*uLl+(1./(1+(sp./sl).^2)).*uPr;

%Sort predictions 
for k=1:F.g1.numL
    for j=1:F.g2.numL
        for i=1:F.g3.numL
            %get mean estimate and its std
            uPo.g1g2g3{i,j,k}=uPo.thisT(F.g1g2g3.posLli{i,j,k});
            uPo.mean(i,j,k)=mean(uPo.g1g2g3{i,j,k});
            uPo.std(i,j,k)=std(uPo.g1g2g3{i,j,k}); 
            %production
            sEs(i,j,k)=uPo.std(i,j,k)+sM;
        end
    end
end

%store the model's information
pred.mean=uPo.mean;
pred.std=sEs;

%.............
%calculate R^2
%.............
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

%...............................
%Calculate std of fit parameters
%...............................
function stdPa=GetStdP(Hessian,data,pred,fitP)
%note: complexe number Std means parameters have negative variance

%Calculate parameters' covariance matrix
%d.covar=inv(jacobian'*jacobian);
d.covar=inv(Hessian); %a matrix

%Calculate noise variance
residual=(data-pred)';
noiseVariance=(residual*residual')/(numel(data) -  numel(fitP));%a scalar

%Calculate std of model parameters
stdPa=diag(sqrt(noiseVariance*d.covar))';

%...........................
%Get prior's representation
%...........................
function Sp=GetPrior1(freeP,factor)

%INPUTS
%freeP: In descending order
%factor: correspond to the freeP

%Convert cell input to matrix
if iscell(freeP)==1
    freeP=[freeP{:}];
end

%Remove cases of infinite prior from the analysis
if freeP(:,1)==+inf 
    freeP(1)=[];
end

%Get the priors' true and estimated values
Sp.exp.raw=factor;
Sp.estimated.raw=freeP;

%Get normalized priors' width (std of weaker prior/ std of prior_i)
for i=1:numel(freeP )
    %Experimental
    Sp.exp.normalized(i)=Sp.exp.raw(1)/Sp.exp.raw(i);%Sw/Ss
    
    %Estimated 
    Sp.estimated.normalized(i)=Sp.estimated.raw(i)/Sp.estimated.raw(1);%Ks/Kw=Sw/Ss
end


%Show summary results
%Variables
fprintf('\n\n\n\n\n %15s %15s \n',...
    'Std.exp.norm',...
    'Std.est.norm')
%Values
for i=1: numel(Sp.exp.normalized)
    fprintf('\n %15i %15f \n',...
        [Sp.exp.normalized(i)';Sp.estimated.normalized(i)']);
end

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
hold all
title ('Subject representations of prior width')
xlabel('{\sigma}_w_e_a_k_e_r _p_r_i_o_r/{\sigma}_p_r_i_o_r','fontsize',12)
ylabel({'K_p_r_i_o_r/K_w_e_a_k_e_r _p_r_i_o_r',...
    'i.e.,{\sigma}_w_e_a_k_e_r _p_r_i_o_r/{\sigma}_p_r_i_o_r'},'fontsize',12)

plot(Sp.exp.normalized,Sp.estimated.normalized,'-ko',...
    'markerfacecolor','k',... 
    'markersize',15,...
    'displayname','subject' );
plot(Sp.exp.normalized,Sp.exp.normalized,'k:',...
    'linewidth',6,'displayname','True ratio');

ylim([1 10])
xlim([0 max(Sp.exp.normalized)+2])
lg=legend('location','Northwest');
box(lg,'off')
axis square
set(gca,'fontsize',12)
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
Sp.estimated.raw=freeP(5:8);

%Show summary results
%Variables
fprintf('\n\n\n\n\n %15s %15s \n',...
    'Std.exp.norm',...
    'Std.est.norm')
%Values
for i=1: numel(Sp.exp.raw)
    fprintf('\n %15i %15f \n',...
        [Sp.exp.raw(i)';Sp.estimated.raw(i)']);
end

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
hold all
title ('Subject representations of prior width')
xlabel('{\sigma}_p_r_i_o_r','fontsize',12)
ylabel('{\sigma}_p_r_i_o_r','fontsize',12)

plot(Sp.exp.raw,Sp.estimated.raw,'-ko',...
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
function Sp=GetPrior5(freeP,stdPa,factor,factor2)

%INPUTS
    %freeP: In descending order
    %factor: factors associated to the freeP

    
%Representation of the prior
%.......................................................---
%Get the priors' true and estimated values
Sp.exp.raw     =factor;
Sp.estimated.raw=freeP(4:7);
Sp.std         =stdPa(4:7);

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
hold all
title ('Subject representations of prior width','fontsize',12)
xlabel('Experimental {\sigma}_p_r_i_o_r','fontsize',12)
ylabel('Estimated {\sigma}_p_r_i_o_r','fontsize',12)
myerrorbar(Sp.exp.raw,Sp.estimated.raw,'yError',Sp.std,...
    'Symbol=o',...
    'Markersize=30',...
    'Color=[0 0 0]');
plot(Sp.exp.raw,Sp.exp.raw,'k:','linewidth',1)
linefit(Sp.exp.raw,Sp.estimated.raw,'k')
set(gca,'xtick',[0:10:90],'xticklabel',[0:10:90],'fontsize',20)

%Representation of the likelihood
%.......................................................---
%Get the llh' true and estimated values
Sl.estimated.raw=freeP(1:3);
Sl.std         =stdPa(1:3);

%Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
figure('color',[1 1 1]);
hold all
title ('Subject representations of llh width','fontsize',12)
xlabel('Coherence','fontsize',12)
ylabel('Estimated {\sigma}_l_l_h','fontsize',12)
myerrorbar(factor2,Sl.estimated.raw,'yError',Sl.std,'Color=[0 0 0]');
xlim([0 1])
ylim([0 Sl.estimated.raw(end)+Sl.std(end)])
set(gca,'xtick',[0.06:0.06:1],'xticklabel',[0.06:0.06:1],'fontsize',20)








%.................
%support functions
%.................
%Calculate the circular statistics of the data
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
% Convert from cartesian coordinates to angles (degree)
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
%Convert from polar to cartesian coordinates
function coord=polar2cartesian(theta,r)
%theta is an angle in degree
%r is the radius of the unit circle
%Coord are in visual angle
%Record angle in degree
theta2.deg=theta;
%Convert from degree to radian
theta2.rad=theta2.deg*pi/180;
%Calculate visual angles coordinates
x=r*cos(theta2.rad);
y=r*sin(theta2.rad);
coord=[x y];
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
%find non missing data
datahere=~isnan(y);
y=y(datahere);

%compute x and y for the linear fit
xfit=x;
P=polyfit(x(datahere),y,1);
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


function axs=makeAxes(F)
%axes width
width=1/(F.g1.numL+1);

%space between axes
gap=(1-(F.g1.numL*width))/(F.g1.numL+1);

%positionne axes
for k=1:F.g1.numL
    axs.position(k,:)=[k*gap+(k-1)*width 0.11 width 0.81];
end




