

%usage

     %(fMRI)

    % v = mrLoadRet([])  
    % v = viewSet(v,'curGroup',2);
    % v = viewSet(v,'curScan',2);
    % o.stimvol = getStimvol(v,'myRandomCoh','taskNum=2','phaseNum=1','segmentNum=2');
    % roi = loadROITSeries(v,'V1.mat');
    % roi = getInstances(v,roi,o.stimvol,'startLag=0','blockLen=2');
    % i = roi{1}.classify.instances;
    % [retval,instancesNew] = slleaveOneOfkmeanInstanceOut(i,2)
 
function [retval,instancesNew] = slleaveOneOfkmeanInstanceOut(instances,k,varargin)

% check arguments
if any(nargin == [0])
  help slkFoldCV
  return
end

% get arguments
type = [];kernelfun = [];kernelargs = [];C=[];fieldName=[];hailString=[];
getArgs(varargin,{'type=fisher','kernelfun=[]','kernelargs=[]','C=[]','fieldName=classify','hailString=[]'});

% see if we are passed in a cell array of rois. If so, then call slkFoldCV
% sequentially on each roi and put the output into the field specified by classField
if isfield(instances{1},fieldName) && isfield(instances{1},'name')
  for iROI = 1:length(instances)
    if ~isfield(instances{iROI}.(fieldName),'instances')
      disp(sprintf('(slkFoldCV) No instances found in %s for %s',fieldName,instances{iROI}.name));
    else
      % put the output into the roi with the field specified by classField
      instances{iROI}.(fieldName).slkFoldCV = slkFoldCV(instances{iROI}.(fieldName).instances,'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,sprintf('hailString=%s%s: ',hailString,instances{iROI}.name));
    end
  end
  retval = instances;
  return
end

% number of classes we are going to classify into
numClasses = length(instances);
% number of repeats we have in each class
for iClass = 1:numClasses
  numReps(iClass) = size(instances{iClass},1);
  numDims(iClass) = size(instances{iClass},2);
  disp(sprintf('(leaveOneOfkmeanInstanceOut) Class %i has %i instances with %i dimensions',iClass,numReps(iClass),numDims(iClass)));
end

% check for dimensions being bad
if length(unique(numDims)) ~= 1
  disp(sprintf('(leaveOneOfkmeanInstanceOut) All instances must have the same number of dimensions'));
  return
end

% save info on how the classification was done
retval.type = type;
retval.classifierParams.kernelfun = kernelfun;
retval.classifierParams.kernelargs = kernelargs;
retval.classifierParams.C = C;

%create k data subsets each containing at least one instance
%of each class. Average instances within each class in each dataset.
%use average instances. Averaging instances by class within test set
%does not break cross-validation because no information is shared 
%between training and test sets. 
%for now we average numReps/k instances contiguous within each class
%but we might want to average over all possible numReps/k instances 
%combinations, and select classification most likely predicted class
%for each test average instance based on vote.

%check that there are enough instances for k fold
if any(numReps/k < 1)
  warning('%s %i %s %i \n','(leaveOneOfkmeanInstanceOut) You need more instances for a',k, '-fold CV. Reduce k to at least',min(numReps))
  return
end

%create k average instance by class
iPerSet = floor(numReps/k); 
for i = 1 : numClasses
  
  sets{i}  = 1 : iPerSet(i) : numReps(i);
  sete{i}  = iPerSet(i) : iPerSet(i) : numReps(i);
  
  %new instances are calculated by 
  %averaging floor(numReps/k) successive instances within
  %the class
  for ki = 1 : k-1
    instancesNew{i}(ki,:) = mean(instances{i}(sets{i}(ki):sete{i}(ki),:),1);    
  end
  
  %the last instance is averaged over all remaining instances
  %within the class
  instancesNew{i}(k,:) = mean(instances{i}(sets{i}(k):end,:),1);
end

% new number of repeats we have in each class
fprintf('%s \n','(leaveOneOfkmeanInstanceOut) Number of repeats after instance averaging:')
for iClass = 1:numClasses
  numReps(iClass) = size(instancesNew{iClass},1);
  numDims(iClass) = size(instancesNew{iClass},2);
  disp(sprintf('(leaveOneOfkmeanInstanceOut) Class %i has %i instances with %i dimensions',iClass,numReps(iClass),numDims(iClass)));
end

%leave-One-new Instance (averaged ones) Out cross-validation
%cycle through class and new instances (averaged), removing a single
%instance, and building the classifier on the remaining
%instances and testing on that single instance.
disppercent(-1/numClasses,sprintf('(leaveOneOfkmeanInstanceOut) %sPerforming k-fold cross-validation for classification from averaged instances with classifier %s',hailString,retval.type));
for iClass = 1:numClasses
  for iRep = 1:k
    % get the test instance
    testInstance = instancesNew{iClass}(iRep,:);
    % cerate the training instances, by removing just the testInstance
    trainingInstances = instancesNew;
    trainingInstances{iClass} = instancesNew{iClass}(setdiff(1:k,iRep),:);
    % now build the classifier
    thisClassifier = buildClassifier(trainingInstances,sprintf('type=%s',type),'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C);
    % and try to classify the instance
    [retval.whichClass{iClass}(iRep) retval.classifierOut{iClass}(iRep)] = classifyInstance(thisClassifier,testInstance);
    % update disppercent
    disppercent((iClass-1)/numClasses,iRep/numReps(iClass));
  end
  % compute how many were correct for this class
  correctByClass(iClass) = sum(retval.whichClass{iClass}==iClass);
  % and compute the confusion matrix row for this class
  for jClass = 1:numClasses
    retval.confusionMatrix(iClass,jClass) = sum(retval.whichClass{iClass}==jClass)/numReps(iClass);
  end
end

% now make into percent correct
retval.correct = sum(correctByClass)/sum(numReps);
disppercent(inf,sprintf('(leaveOneOfkmeanInstanceOut) %s%s classifier produced %0.2f%% correct and',hailString,retval.type,retval.correct*100));
fprintf('\n')
retval.correctSTE = sqrt(retval.correct*(1-retval.correct)/sum(numReps));



