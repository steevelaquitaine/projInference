

%slfmriSimInstances.m
%
%
% author: steeve laquitaine
%   date: 151202
%purpose: simulate instances to test classification
%         c{1}.classify.instances{1} are instances for class 1 roi 1
%         c{1}.classify.instances{2} are instances for class 2 roi 1
%         etc,...
%
%  usage:
%
% ex: 1
%           c=slfmriSimInstances({'simV1','simMT'},2,[20 30],[10 20])
%           c=leaveOneOut(c,'permutation=1');
%
%
% ex:2
%           %simulate two 2D-Gaussian-clusters of instances
%           c = slfmriSimInstances({'simV1'},6,[],[50 50 50 50 50 50],'type=a2Dclusters')                   
%           %c = c{1}.classify.instances;
%           %c = leaveOneOut(c);
%           %plot
%           classs = linspecer(6);
%           for i = 1 : length(c{1}.classify.instances)
%              hold on;
%              plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),'.','color',classs(i,:))
%           end
% 
% ex:3
%           %2D-Gaussian-clusters 
%           c=slfmriSimInstances({'simV1','simMT'},2,[],[10 15],'type=a2Dclusters')         
%           c=leaveOneOut(c,'balancByBootSt=1');           
%           %plot
%           classs={'r','b'};
%           nClass = length(c{1}.classify.instances);
%           for i=1:nClass
%              hold on
%              plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),...
%              '.','color',classs{i})
%           end
%
% ex4:
%           %2D-Gaussian-clusters 
%           c=slfmriSimInstances({'simV1','simMT'},2,[],[10 15],'type=a2Dclusters')
%           %classify
%           c=leaveOneOut(c,'balancByRemovI=1');
%           %plot
%           classs={'r','b'};
%           nClass=length(c{1}.classify.instances);
%           for i=1:nClass
%               hold on
%               plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),...
%               '.','color',classs{i})
%           end
%
% ex5:
%           %2D-Gaussian-clusters         
%           c=slfmriSimInstances({'simV1','simMT'},2,[],[20 20],'type=a2Dclusters')
%           %classify
%           c=kFold(c,'numFolds=10');
%           %plot
%           classs={'r','b'};
%           nClass=length(c{1}.classify.instances);
%           for i=1:nClass
%               hold on
%               plot(c{1}.classify.instances{i}(:,1),c{1}.classify.instances{i}(:,2),...
%               '.','color',classs{i})
%           end

function c = slfmriSimInstances(rois,nClasses,nvoxs,numInst,varargin)

type=[];
getArgs(varargin,{'type=rnd'});

%# of rois
nRois = length(rois);

%(default) case random instance values
if strcmp(type,'rnd')
    %make instances for each roi and
    %class
    for roi = 1 : nRois
        %name roi
        c{roi}.name = rois{roi};
        for class = 1 : nClasses
            c{roi}.classify.instances{class} = rand(numInst(class),nvoxs(roi));
        end
    end
end

%case 2D cluster of n classes
%2 voxels
mu = 1:nClasses;
sigma = eye(nClasses,nClasses);
%sigma = [10 4; 1 5];
if strcmp(type,'a2Dclusters')
    %make instances for each roi and
    %class
    for roi = 1 : nRois
        %name roi
        c{roi}.name = rois{roi};
        for class = 1 : nClasses
            m = mu + (class - 1)*5;     
            c{roi}.classify.instances{class} = mvnrnd(m,sigma,numInst(class))*100;
        end
    end
end