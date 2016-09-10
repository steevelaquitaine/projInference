
 %author: steeve laquitaine
   %date: 140410
  %usage:
%purpose: getInstancesVarnames


function [icat,icatnames]=getInstancesVarnames(i,stimNames)
Ibkp=[];
icat=[];
for j=1:length(i)
    
    %variables name
    inames=stimNames{j};
    inames=cellstr(inames(ones(size(i{j},1)),:));
    icatnames=[Ibkp;inames];
    
    %concatenate instances, icat (n trials,m voxels)
    icat=[icat;i{j}];
end
    