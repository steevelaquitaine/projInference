

%slfMRIgetAnalyses.m
%
%purpose: store the list of analysis parameters that will be run

%get analysis
function o = slfMRIgetAnalyses(o,varargin)
%get arg
varargin = varargin{:};
%calculate accuracy each time step
if sum(strcmp(varargin,'accuracyByTime'))
    o.myAnalysis = {'accByTime'};
elseif sum(strcmp(varargin,'accuracyAtTime'))
    %accuracy at specific time
    o.myAnalysis = {'accAtTime'};
    pos = find(strcmp(varargin,'accuracyAtTime')) + 1;
    o.volsToClass = varargin{pos};                        %[start end] vols to average for classification
end
%if testing 
if sum(strcmp(varargin,'test'))
    o.test=1;
else
    o.test=0;
end
%--------- cross-validation -----------
%leave-one instance out
if sum(strcmp(varargin,'leaveOneOut'))
    o.myAnalysis{end+1} = 'leaveOneOut';
    o.leaveOneOut = 1;
    o.crossVal = 'leaveOneOut';
else
    o.leaveOneOut = 0;
end
%k-Fold
if sum(strcmp(varargin,'kFoldCV'))
    o.myAnalysis{end+1} = 'kFoldCV';
    o.kFoldCV = 1;
    o.crossVal = 'kFoldCV';
else
    o.kFoldCV = 0;
end
%leave-one k mean instance out
%group instances into K average instances by class and leave-one-out
%classification. The idea is to average out past trial Bold confound.
if sum(strcmp(varargin,'slleaveOneOfkmeanInstanceOut'))
    o.myAnalysis{end+1} = 'slleaveOneOfkmeanInstanceOut';
    pos = find(strcmp(varargin,'slleaveOneOfkmeanInstanceOut')) + 1;
    o.kInstances = varargin{pos};
    o.crossVal = 'slleaveOneOfkmeanInstanceOut';
end
%--------------classifier ------------
%check first
if ~any(strcmp(varargin,{'fisher'})) && ~any(strcmp(varargin,{'svm'}))
    fprintf('%s \n','(slgetAnalysis) Please input classifier: "fisher" or "svm"')
    dbstack
    keyboard
end
if sum(strcmp(varargin,'svm'))
    o.myAnalysis{end+1} = 'svm';
    o.type = 'svm';
end
if sum(strcmp(varargin,'fisher'))
    o.myAnalysis{end+1} = 'fisher';
    o.type = 'fisher';
end
%---------- denoising -------------------
if sum(strcmp(varargin,'r2cutoff'))
    o.myAnalysis{end+1} = 'r2cutoff';
    pos = find(strcmp(varargin,'r2cutoff')) + 1;
    o.r2thres = varargin{pos};    
end
%regress out a variable (not yet correct)
if any(strcmp(varargin,'regOutVar2'))
    o.myAnalysis{end+1} = 'regOutVar2';
    pos = find(strcmp(varargin,'regOutVar2')) + 1;
    o.regOutVar2 = varargin{pos};
end

%----------balancing dataset -------------
if any(strcmp(varargin,'balancByRemovI=1'))
    o.myAnalysis{end+1} = 'balancByRemovI';    
    o.dataBalancing = 'balancByRemovI=1';
elseif any(strcmp(varargin,'balancByRemovI=0'))
    o.dataBalancing = 'balancByRemovI=0';
end
if any(strcmp(varargin,'balancByBootSt=1'))
    o.myAnalysis{end+1} = 'balancByBootSt';    
    o.dataBalancing = 'balancByBootSt=1';
elseif any(strcmp(varargin,'balancByBootSt=0'))
    o.dataBalancing = 'balancByBootSt=0';
end

%---------- chance level ----------------
%calculate chance by permutation
if sum(strcmp(varargin,'chanceByPerm'))
    o.myAnalysis{end+1} = '_chanceByPerm';
    o.chanceByPerm = 1;
else
    o.chanceByPerm = 0;
end
%balancing
if any(strcmp(varargin,'permutationUnBal=1'))
    o.myAnalysis{end+1} = 'permutationUnBal';    
    o.chanceBalancing = 'permutationUnBal=1';
end
if any(strcmp(varargin,'permutationBal=1'))
    o.myAnalysis{end+1} = 'permutationBal';    
    o.chanceBalancing = 'permutationBal=1';
end

%store varargin
o.myvarg = varargin;