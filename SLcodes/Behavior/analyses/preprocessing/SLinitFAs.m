
%SLinitFAs.m
%
%       $Id: SLinitFAs.m 750 2012-09-08 07:12:46Z steeve $
%     usage: SLinitFAs(databank,FAs)
%        by: steeve laquitaine
%      date: 140531 21:34:00

%   purpose: find the factors in a data matrix
%          

%%Set the factors (predictors) of the analysis
function [FA,est_dir]=SLinitFAs(databank,FAs)

%get data
%--------
%Get est data (e.g., estimated dirs)
est_data   =cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));

%Set factors
%FA 1
FA.g1.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{1}))==1));
FA.g1.nm       =FAs{1};
%FA 2
FA.g2.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{2}))==1));
FA.g2.nm       =FAs{2};
%FA 3
FA.g3.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{3}))==1));
FA.g3.nm       =FAs{3};


%CASE 1: Prior is unimodal
%-------------------------
%Get and order levels of each group
%group 1 (e.g., priors)
FA.g1.lvlsnm=unique(FA.g1.thisT); %names
FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend'); %order
FA.g1.lvlsnb=numel(FA.g1.lvlsnm);   %number
clear i
for i=1:FA.g1.lvlsnb
    index.g1.lvli(i)={find(FA.g1.thisT==FA.g1.lvlsnm(i))};
end
%group 2 (e.g., coherences)
FA.g2.lvlsnm=unique(FA.g2.thisT); %names
FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend'); %order
FA.g2.lvlsnb=numel(FA.g2.lvlsnm);   %number
for i=1:FA.g2.lvlsnb
    index.g2.lvli(i)={find(FA.g2.thisT==FA.g2.lvlsnm(i))};
end
%group 3 (e.g., displayed dirs)
FA.g3.lvlsnm=unique(FA.g3.thisT); %names
FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'ascend'); %order
FA.g3.lvlsnb=numel(FA.g3.lvlsnm);   %number
for i=1:FA.g3.lvlsnb
    index.g3.lvli(i)={find(FA.g3.thisT==FA.g3.lvlsnm(i))};
end


%CASE 2, Prior is bimodal
%------------------------
%factor 1 are the prior modes sorted based on how far they are from each
%other: the distance between the modes.

%if factor 1 is the prior
clear i
if strcmp(FA.g1.nm,'priormodes')

    %record the distance between the modes of each bimodal prior
    FA.g1.distance_mode=FA.g1.thisT(:,2)-FA.g1.thisT(:,1);
    FA.g1.lvlsnm=unique(FA.g1.distance_mode); 
    FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend');
    FA.g1.lvlsnb=numel(FA.g1.lvlsnm);

    %find trials position of each prior condition
    for i=1:FA.g1.lvlsnb
        index.g1.lvli(i) = {find(FA.g1.distance_mode==FA.g1.lvlsnm(i))};
    end
end

%if factor 2 is the prior
clear i
if strcmp(FA.g2.nm,'priormodes')
    
    %record the distance between the modes of each bimodal prior
    FA.g2.distance_mode=FA.g2.thisT(:,2)-FA.g2.thisT(:,1);
    FA.g2.lvlsnm=unique(FA.g2.distance_mode);
    FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend');
    FA.g2.lvlsnb=numel(FA.g2.lvlsnm);
    
    %find trials position of each prior condition
    for i=1:FA.g2.lvlsnb
        index.g2.lvli(i)={find(FA.g2.distance_mode==FA.g2.lvlsnm(i))};
    end
    
end

%if factor 3 is the prior
clear i
if strcmp(FA.g3.nm,'priormodes')

    %record the distance between the modes of each bimodal prior
    FA.g3.distance_mode=FA.g3.thisT(:,2)-FA.g3.thisT(:,1);
    FA.g3.lvlsnm=unique(FA.g3.distance_mode); 
    FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'descend');
    FA.g3.lvlsnb=numel(FA.g3.lvlsnm);

    %find trials position of each prior condition
    for i=1:FA.g3.lvlsnb
        index.g3.lvli(i)={find(FA.g3.distance_mode==FA.g3.lvlsnm(i))};
    end
end

%Calculate coordinates of average estimated directions for each condition
%organize data in following order: group 1(subplots) - group 2(colors) -
%group 3(x-axis). Each cell contains repetitions of a condition.
for k=1:FA.g1.lvlsnb
    for j=1:FA.g2.lvlsnb
        for i=1:FA.g3.lvlsnb
            index.g1g2g3(j,i,k)={intersect( intersect( index.g1.lvli{k},index.g2.lvli{j} ),index.g3.lvli{i} )};
            
            %calculate statistics of estimated dirs (j:group 2, i:group 3, k:group 1)
            est_dir{j,i,k}=SLstatcircular(est_data(index.g1g2g3{j,i,k},:));
            
        end
    end
end
