
%slfmriClassify2.m
%
%    Author: steeve laquitaine
%      date: 151225
%   purpose: decoding various features with various options from ROI 
%            voxels Bold pattern
%
%inputs:
%
%       'CalcInstances'(or 'loadSavedInstances')
%       'slReSampInst': re-sample instances when few only
%                       set random values for missing classes. Do not run
%                       when classes are missing !!
%               'test': just testing: save  in a test folder
%
%Loading:
%         'loadSavedROI' : load existing roi TSeries file loaded with
%                          loadROITSeries which is time expensive
%
%cross-validation:
%          'leaveOneOut' : train on all instance-1 and test on left over
%              'kFoldCV' : divide into 10 folds train on 9 folds test on 1
%                          left over
%
%
%Chance level:
%        'chanceByPerm' : calculate chance label from N-permuted dataset
%            'nPerm=100': to get chance level form 100 dataset of permuted
%                         classes
%     'permutationUnBal': permute instances classes
%       'permutationBal': equalize number of instances and permute classes
%
%denoising:
%       'r2cutoff',0.03: get rid of voxels with low r2 in an event-related
%                        analysis. You have to indicate before the name of 
%                        the er file that should be loaded (e.g., o.erAnalFile =
%                        'erAnalMotionByAll.mat')
%
%outputs:
%--------



function slfmriClassify2(o,myClass,varargin)

t0=tic;

%get analysis parameters
o = slfmriGetClassifAnalysisParams(o,myClass,varargin{:});

%stack instances over sessions
[o,si,s] = slStackSessionInstances(o,varargin{:});

%classify by roi
nVar = length(o.taskCond);
for roii = 1 : o.nRois
    fprintf('%s \n',['(leaveOneOut)' o.myROIname{roii}] )
    
    %path
    cd([o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]])
    
    %init   
    iStacked  = s.iStacked{roii};
    izStacked = s.izStacked{roii};

    %---------------------------- for debugging ----------------
    %simulate instances for test
    %c=simInstances({'lV1'},nVar,[],repmat(25,nVar,1),'type=a2Dclusters');
    %iStacked=repmat({c{1}.classify.instances},o.nRois,1);
    %izStacked=preprocessInstances(izStacked{1},'zscore=1');
    %izStacked=repmat({izStacked},o.nRois,1);
    %---------------------------- for debugging ----------------
    
    %---------------------------- for debugging ----------------
    %increase small dataset by resampling
    %size of dataset.just for test: might break
    %the training/test isolation requirement.
    %also I set to rand the classes when instances
    %that have no instances
    %this produces high %correct
    %if slIsInput(varargin,'slReSampInst')
    %    [iStacked,izStacked] = slReSampInst(iStacked,izStacked);
    %end
    %---------------------------- for debugging ----------------            
    
    fprintf('\n %s \n \n','------ Classification (actual data) --------')
    c.myClasf.raw.fullSg = [];
    c.myClasf.Zsc.fullSg = [];
    
    %checking when empty classes
    for ij = 1: length(iStacked)
        emptyClass(ij) = isempty(iStacked{ij});
    end
    iStacked2 = iStacked;
    izStacked2 = izStacked;
    iStacked2(emptyClass) = [];
    izStacked2(emptyClass) = [];
          
    %case leaveoneout
    if slIsInput(varargin,'leaveOneOut')        
        
        c.myClasf.raw.fullSg        = leaveOneOut(iStacked2,['type=' o.type],o.dataBalancing);%raw                
        c.myClasf.raw.fullSgClass   = s.variable(~emptyClass);
        c.myClasf.raw.TheoretChance = 1 / length(s.variable(~emptyClass));
        %c.myClasf.Zsc.fullSg = leaveOneOut(izStacked2,['type=' o.type],o.dataBalancing);%zs
        c.myClasf.Zsc.fullSgClass   = s.variable(~emptyClass);
        c.myClasf.Zsc.TheoretChance = 1 / length(s.variable(~emptyClass));
    end
    
    %case k-fold
    if slIsInput(varargin,'kFoldCV')
        fprintf('%s \n','kFoldCV')
        fprintf('%s \n \n','-----------')        
        c.myClasf.raw.fullSg = kFold(iStacked2,'stratified','numFolds=10',['type=' o.type],o.dataBalancing);%raw
        disp(sprintf('(kFold) %s classifier produced %0.2f%% correct',o.type,c.myClasf.raw.fullSg.correct*100));
        
        %c.myClasf.Zsc.fullSg = kFold(izStacked2,'stratified','numFolds=10',['type=' o.type],o.dataBalancing);%zs
        %disp(sprintf('(kFold) %s classifier produced %0.2f%% correct',o.type,c.myClasf.Zsc.fullSg.correct*100));
        
        c.myClasf.raw.fullSgClass = s.variable(~emptyClass);
        c.myClasf.raw.TheoretChance = 1 / length(s.variable(~emptyClass));        
        c.myClasf.Zsc.fullSgClass = s.variable(~emptyClass);
        c.myClasf.Zsc.TheoretChance = 1 / length(s.variable(~emptyClass));
    end 
    fprintf('\n')
    fprintf('%s \n \n','---- (Classification at chance) permuted data ---')
    %classify permuted data N times (typically N=100)
    %36 min (50 perms,1 ROI) (Markus Ojala, see Permutation 
    %test for stuyding classifer performance, Journal of 
    %Machine learning, 2010) 
    o.instances    = iStacked2;
    o.Zscinstances = izStacked2;
    [c,o] = slGetClassifChance(o,c);

    %fprintf('%s \n','------ Classification (balanced data) ---------')
    %uncomment. accuracy doesn't change much when we balance the dataset
    %balance data (removing) and classify
    c.myClasf.raw.fullSgShufBal = [];
    c.myClasf.Zsc.fullSgShufBal = [];
    %c.myClasf.raw.fullSgShufBal = leaveOneOut(iStacked{roii},'balancByRemovI=1');%raw
    %c.myClasf.Zsc.fullSgShufBal = leaveOneOut(izStacked{roii},'balancByRemovI=1');%zscore
    
    %save results classification
    [o,c,s,roii] = slSaveClassiRes(o,c,s,roii);
    
    %duration
    duration=toc(t0);
    fprintf('%s %.2f \n','One more roi --------------',duration)
end

%Go to analysis path
% o.analysisPath = [o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]];
% cd(o.analysisPath)

% 
% %--------------- nested -------------
% %get analysis
% function o = slgetAnalysis(o,varargin)
% %get arg
% while size(varargin)==1
%     varargin = varargin{:};
% end
% %calculate accuracy each time step
% if sum(strcmp(varargin,'accuracyByTime'))
%     o.myAnalysis = {'accByTime'};
% elseif sum(strcmp(varargin,'accuracyAtTime'))
%     %accuracy at specific time
%     o.myAnalysis = {'accAtTime'};
%     pos = find(strcmp(varargin,'accuracyAtTime')) + 1;
%     o.volsToClass = varargin{pos};                        %[start end] vols to average for classification
% end
% %if testing 
% if sum(strcmp(varargin,'test'))
%     o.test=1;
% else
%     o.test=0;
% end
% %--------- cross-validation -----------
% %leave-one instance out
% if sum(strcmp(varargin,'leaveOneOut'))
%     o.myAnalysis{end+1} = 'leaveOneOut';
%     o.leaveOneOut = 1;
%     o.crossVal = 'leaveOneOut';
% else
%     o.leaveOneOut = 0;
% end
% %k-Fold
% if sum(strcmp(varargin,'kFoldCV'))
%     o.myAnalysis{end+1} = 'kFoldCV';
%     o.kFoldCV = 1;
%     o.crossVal = 'kFoldCV';
% else
%     o.kFoldCV = 0;
% end
% %leave-one k mean instance out
% %group instances into K average instances by class and leave-one-out
% %classification. The idea is to average out past trial Bold confound.
% if sum(strcmp(varargin,'slleaveOneOfkmeanInstanceOut'))
%     o.myAnalysis{end+1} = 'slleaveOneOfkmeanInstanceOut';
%     pos = find(strcmp(varargin,'slleaveOneOfkmeanInstanceOut')) + 1;
%     o.kInstances = varargin{pos};
%     o.crossVal = 'slleaveOneOfkmeanInstanceOut';
% end
% %--------------classifier ------------
% %check first
% if ~any(strcmp(varargin,{'fisher'})) && ~any(strcmp(varargin,{'svm'}))
%     fprintf('%s \n','(slgetAnalysis) Please input classifier: "fisher" or "svm"')
%     dbstack
%     keyboard
% end
% if sum(strcmp(varargin,'svm'))
%     o.myAnalysis{end+1} = 'svm';
%     o.type = 'svm';
% end
% if sum(strcmp(varargin,'fisher'))
%     o.myAnalysis{end+1} = 'fisher';
%     o.type = 'fisher';
% end
% %---------- denoising -------------------
% if sum(strcmp(varargin,'r2cutoff'))
%     o.myAnalysis{end+1} = 'r2cutoff';
%     pos = find(strcmp(varargin,'r2cutoff')) + 1;
%     o.r2thres = varargin{pos};    
% end
% %regress out a variable (not yet correct)
% if any(strcmp(varargin,'regOutVar2'))
%     o.myAnalysis{end+1} = 'regOutVar2';
%     pos = find(strcmp(varargin,'regOutVar2')) + 1;
%     o.regOutVar2 = varargin{pos};
% end
% 
% %----------balancing dataset -------------
% if any(strcmp(varargin,'balancByRemovI=1'))
%     o.myAnalysis{end+1} = 'balancByRemovI';    
%     o.dataBalancing = 'balancByRemovI=1';
% elseif any(strcmp(varargin,'balancByRemovI=0'))
%     o.dataBalancing = 'balancByRemovI=0';
% end
% if any(strcmp(varargin,'balancByBootSt=1'))
%     o.myAnalysis{end+1} = 'balancByBootSt';    
%     o.dataBalancing = 'balancByBootSt=1';
% elseif any(strcmp(varargin,'balancByBootSt=0'))
%     o.dataBalancing = 'balancByBootSt=0';
% end
% 
% %---------- chance level ----------------
% %calculate chance by permutation
% if sum(strcmp(varargin,'chanceByPerm'))
%     o.myAnalysis{end+1} = '_chanceByPerm';
%     o.chanceByPerm = 1;
% else
%     o.chanceByPerm = 0;
% end
% %balancing
% if any(strcmp(varargin,'permutationUnBal=1'))
%     o.myAnalysis{end+1} = 'permutationUnBal';    
%     o.chanceBalancing = 'permutationUnBal=1';
% end
% if any(strcmp(varargin,'permutationBal=1'))
%     o.myAnalysis{end+1} = 'permutationBal';    
%     o.chanceBalancing = 'permutationBal=1';
% end
% 
% %store varargin
% o.myvarg = varargin;



%resample instances to increase
%sample size
function [iStacked,izStacked] = slReSampInst(iStacked,izStacked)

%get # vox
nClasses = length(iStacked{:});
nvox=[];
for cl = 1 : nClasses
    nvox = max([nvox size(iStacked{:}{cl})]);
end

nReSamp = 33;
for cl = 1 : nClasses
    nI = size(iStacked{:}{cl},1);
    %if no data
    if isempty(iStacked{:}{cl})
        iStacked{:}{cl} = rand(20,nvox);%raw
        izStacked{:}{cl} = rand(20,nvox);%zscore
    else
        nReSampix = randsample(1:nI,nReSamp,'true');
        iStacked{:}{cl} = iStacked{:}{cl}(nReSampix,:);%raw
        izStacked{:}{cl} = izStacked{:}{cl}(nReSampix,:);%zscore
    end
end

%save classification data
function [o,c,s,roii] = slSaveClassiRes(o,c,s,roii)

%convert cell class to string
o.savedClass = o.myClass;
while iscell(o.savedClass)
    o.savedClass = [o.savedClass{:};];
end
%if not testsing save data in structured dir
if o.test~=1
    %backup
    if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
    cd(o.myROIname{roii})
    o.myAnalPath2{roii} = pwd;
    %archive existing file
    list = dir('*.mat');
    if ~isempty(list)
        if isempty(dir('Archive')); mkdir Archive; end
        for i=1:length(list); movefile(list(i).name,'Archive/'); end
    end
    save(['ClassifStckSess' o.savedClass '_' o.myROIname{roii} '.mat'],'o','c','s')
    %save chance accuracy results
    %in a different file
    chance.raw = c.myClasf.raw.fullSgShuf;
    chance.Zsc = c.myClasf.Zsc.fullSgShuf;
    save(['ClassifStckSessChance' o.savedClass '_' o.myROIname{roii} '.mat'],'chance')
    clear chance;
    cd ..

elseif o.test==1
    %otherwise save in test folder
    %backup
    if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
    cd(o.myROIname{roii})
    mkdir test; cd test    
    o.myAnalPath2{roii} = pwd;
    %archive existing file
    list = dir('*.mat');
    if ~isempty(list)
        if isempty(dir('Archive')); mkdir Archive; end
        for i=1:length(list); movefile(list(i).name,'Archive/'); end
    end
    save(['ClassifStckSess' o.savedClass '_' o.myROIname{roii} '.mat'],'o','c','s')
    %save chance accuracy results
    %in a different file
    chance.raw = c.myClasf.raw.fullSgShuf;
    chance.Zsc = c.myClasf.Zsc.fullSgShuf;   
    save(['ClassifStckSessChance' o.savedClass '_' o.myROIname{roii} '.mat'],'chance')
    clear chance;
end

%cleanup
if isfield(c.myClasf.raw,'fullSgShuf')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShuf');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShuf');
end
if isfield(c.myClasf.raw,'fullSgShufBal')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShufBal');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShufBal');
end

%calculate chance level
function [c,o] = slGetClassifChance(o,c)

%inputs
%   o.crossVal (e.g., 'leaveOneOut' or 'kFold')
%   o.nPerm    (e.g., 2)

%init
c.myClasf.raw.fullSgShuf = [];
c.myClasf.Zsc.fullSgShuf = [];

if o.chanceByPerm==1
    %case leaveOneOut
    if strcmp(o.crossVal,'leaveOneOut')
        if ~isempty(o.nPerm)
            fprintf('%s %i %s \n','(slfmriClassify) Running ',o.nPerm,'permutations...')
            nPerm2 = ['nPerm=' num2str(o.nPerm)];
            type = ['type=' o.type];
            balancing = o.chanceBalancing;
            c.myClasf.raw.fullSgShuf = leaveOneOutNpermut(o.instances,nPerm2,type,balancing);
            %c.myClasf.Zsc.fullSgShuf = leaveOneOutNpermut(o.Zscinstances,nPerm2,type,balancing);
        end
    end
    %case ~10-Fold
    if strcmp(o.crossVal,'kFoldCV')                
        if ~isempty(o.nPerm)
            fprintf('%s %i %s \n','(slfmriClassify) Running',o.nPerm,'permutations to get chance.')
            nPerm2 = ['nPerm=' num2str(o.nPerm)];
            type = ['type=' o.type];
            balancing = o.chanceBalancing;
            c.myClasf.raw.fullSgShuf = kFoldNpermut(o.instances,'stratified',nPerm2,'numFolds=10',type,balancing);
            %c.myClasf.Zsc.fullSgShuf = kFoldNpermut(o.Zscinstances,nPerm2,'numFolds=10',type,balancing);
        end
    end
end



