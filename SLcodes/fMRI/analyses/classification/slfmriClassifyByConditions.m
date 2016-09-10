%slfmriClassifyByConditions.m
%
%
% author: steeve laquitaine
%purpose: classify motion directions from the brain activity for data
%         sorted by a condition (e.g., "switching") to determine whether the 
%         sensory evidence is represented or not when the subjects switch 
%         to his prior.
%         The # of instances used for classification is equalized
%         between the conditions to allow a fair comparison of the
%         classification accuracy between the conditions.
%
%  usage :
% 
%           %setup
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%           params = {'accuracyAtTime',[7 14],'loadSavedROI','CalcInstances',...
%                     'leaveOneOut','fisher','balancByRemovI=1',...
%                     'numInstances',8,'loadInstancesByCond'};           
%           myConds = {{'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06'},....
%                      {'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06'}};
% 
%           %bootstrap and classify
%           parfor i = 1 : 30
%               [acc,accSTE,cbyc{i},sbyc{i},os] = slfmriClassifyByConditions(o,myConds,params);
%                accAll(i,:) = cell2mat(acc);
%                accSTEAll(i,:) = cell2mat(accSTE);
%           end
% 
%           %Get accuracy stats
%           [err_margin1,ci1,m1,s1] = slMakeCI(accAll(:,1),0.95)
%           [err_margin2,ci2,m2,s2] = slMakeCI(accAll(:,2),0.95)
% 
%
%inputs
%
%         'randomseed': if 'randomseed' if input make sure that a seed
%                       variable s.mat is present in the current directory.
%                       This seed has been created with s = rng and saved.
%                       Otherwise if absent just sample randomly instances
%                       each time it is run.
%
%
%The code sometimes generate random values. A random seed is used and can
%be changed.

function [acc,accSTE,cbyc,sbyc,o,obyc] = slfmriClassifyByConditions(o,myConds,params)

paramstmp = params;

t0 = tic ;
%Stack instances and save instance database
if any(strcmp(params,'CalcInstancesByCond')==1)    
    [obyc,sbyc,o] = slfmriStackAndSaveSessionInstances(o,myConds,params);    
    %or load
elseif any(strcmp(params,'loadInstancesByCond')==1)    
    cd([o.classifResultPath '/slfmriClassifyByConditions/' o.myROIname{1} '/' cell2mat([myConds{:};]) '/'])
    load('dataInstancesByCond.mat')
end

%set equal # of instances per condition
if any(strcmp(paramstmp,'randomseed')); 
    o.randomStatus = 'randomseed';  params{end+1} = 'randomseed';
else ; o.randomStatus = [];  end   
for myc = 1 : length(myConds) 
    sbyc{myc}.iStacked{1} = slsetNumInstancesByClass(sbyc{myc}.iStacked{1},o.numInstances,o.randomStatus);
end

%classify for each condition for each roi
for myc = 1 : length(myConds)    
          
    %classify by roi
    for roii = 1 : obyc{myc}.nRois
                
        %info
        fprintf('%s \n',['Classification of data from ' obyc{myc}.myROIname{roii}] )        
        
        %move to saved data path
        cd([o.stckPath 'slStckAnalyses/' obyc{myc}.myGroup '/classif/' obyc{myc}.savedClass '/' [obyc{myc}.myAnalysis{:}]])
        
        %get roi session-stacked instances
        iStacked  = sbyc{myc}.iStacked{roii};
        izStacked = sbyc{myc}.izStacked{roii};
        c.myClasf.raw.fullSg = [];
        c.myClasf.Zsc.fullSg = [];
        
        %checking empty classes
        for ij = 1 : length(iStacked)
            emptyClass(ij) = isempty(iStacked{ij});
        end
        
        iStacked2 = iStacked;
        izStacked2 = izStacked;
        iStacked2(emptyClass) = [];
        izStacked2(emptyClass) = [];
        
        %case leave-one-out
        if slIsInput(params,'leaveOneOut')
            c.myClasf.raw.fullSg        = leaveOneOut(iStacked2,['type=' obyc{myc}.type],obyc{myc}.dataBalancing);
            c.myClasf.raw.fullSgClass   = sbyc{myc}.variable(~emptyClass);
            c.myClasf.raw.TheoretChance = 1 / length(sbyc{myc}.variable(~emptyClass));
            %c.myClasf.Zsc.fullSg = leaveOneOut(izStacked2,['type=' o.type],o.dataBalancing);%zs
            c.myClasf.Zsc.fullSgClass   = sbyc{myc}.variable(~emptyClass);
            c.myClasf.Zsc.TheoretChance = 1 / length(sbyc{myc}.variable(~emptyClass));
        end
        
        %case k-fold
        if slIsInput(params,'kFoldCV')
            fprintf('%s \n','kFoldCV')
            fprintf('%s \n \n','-----------')
            c.myClasf.raw.fullSg = kFold(iStacked2,'stratified','numFolds=10',['type=' obyc{myc}.type],obyc{myc}.dataBalancing);%raw
            disp(sprintf('(kFold) %s classifier produced %0.2f%% correct',obyc{myc}.type,c.myClasf.raw.fullSg.correct*100));
            
            %c.myClasf.Zsc.fullSg = kFold(izStacked2,'stratified','numFolds=10',['type=' o.type],o.dataBalancing);%zs
            %disp(sprintf('(kFold) %s classifier produced %0.2f%% correct',o.type,c.myClasf.Zsc.fullSg.correct*100));            
            c.myClasf.raw.fullSgClass = sbyc{myc}.variable(~emptyClass);
            c.myClasf.raw.TheoretChance = 1 / length(sbyc{myc}.variable(~emptyClass));
            c.myClasf.Zsc.fullSgClass = sbyc{myc}.variable(~emptyClass);
            c.myClasf.Zsc.TheoretChance = 1 / length(sbyc{myc}.variable(~emptyClass));
        end
        fprintf('\n')
        fprintf('%s \n \n','---- (Classification at chance) permuted data ---')
        %classify permuted data N times (typically N=100)
        %36 min (50 perms,1 ROI) (Markus Ojala, see Permutation
        %test for stuyding classifer performance, Journal of
        %Machine learning, 2010)
        obyc{myc}.instances    = iStacked2;
        obyc{myc}.Zscinstances = izStacked2;
        [c,obyc{myc}] = slGetClassifChance(obyc{myc},c);
        
        %fprintf('%s \n','------ Classification (balanced data) ---------')
        %uncomment. accuracy doesn't change much when we balance the dataset
        %balance data (removing) and classify
        c.myClasf.raw.fullSgShufBal = [];
        c.myClasf.Zsc.fullSgShufBal = [];
        %c.myclasf.raw.fullSgShufBal = leaveOneOut(iStacked{roii},'balancByRemovI=1');%raw
        %c.myclasf.Zsc.fullSgShufBal = leaveOneOut(izStacked{roii},'balancByRemovI=1');%zscore
        
        %save results classification (doesn't work with parfor)
        %[obyc{myc},c,sbyc{myc},roii] = slSaveClassiRes(obyc{myc},c,sbyc{myc},roii);                
        cbyc{myc} = c;
        
        %duration
        duration = toc(t0);
        fprintf('%s %.2f \n','One more roi --------------',duration)
        
        %store accuracies    
        acc{myc} = cbyc{myc}.myClasf.raw.fullSg.correct;
        accSTE{myc} = cbyc{myc}.myClasf.raw.fullSg.correctSTE;
    end
end

fprintf('%s \n','----------------------------------')
for myc = 1 : length(myConds)
    fprintf('%s %i %s %.02f %s %.02f \n','%correct condition: ',myc,': ',cbyc{myc}.myClasf.raw.fullSg.correct, '; s:',cbyc{myc}.myClasf.raw.fullSg.correctSTE)
end
fprintf('%s \n','-------------------------------------')
