
%slGetClassifChance.m
%
%author : steeve laquitaine
%purpose:
%
%  usage :
%
%       [c,o] = slGetClassifChance(o,c)

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