
%slfmriGraphPermutedVoxTunings.m
%
%
% author : steeve laquitaine
%purpose : visualize individual tunings produced by permuting voxel 
%          responses y to stimulus feature x 
%
% usage : 
%
%       [ymean_by_xu,xu] = slfmriGraphPermutedVoxTunings(y,x,100)
%
%inputs :
%       
%       y : voxel responses instances
%       x : stimulus feature (e.g., motion direction)
%   nperm : number of permutations


function [ymean_by_xu,xu,datadegmean,tunStrngth] = slfmriGraphPermutedVoxTunings(y,x,nperm)

num = length(y);

%plot observed tuning
figure('color','w')
[obsymean_by_xu,obsxu] = slCalculateMeanSortedBy1Var(y,x);    
hold on; h = plot(obsxu,obsymean_by_xu,'-ro');
title('Observed tuning')

figure('color','w')
h = []; hh = [];
for i = 1 : nperm       
    
    %permute vox responses
    j = randperm(num);
    yperm = y(j);
    
    %get mean y responses for each x
    %Averaging removes the effect of the non-uniform
    %distribution of x on the voxel feature vector (0.001)
    [ymean_by_xu{i},xu] = slCalculateMeanSortedBy1Var(yperm,x);    
    hold on; h(i) = plot(xu,ymean_by_xu{i},'-ro');
    if i>1
        set(h(1:end-1),'color','b'); set(h(i),'color','r')
    end
    drawnow 
            
    %von mises fit    
    [datadegmean(i),tunStrngth(i)] = slfmriGetVoxTuningVMfitGridSearch(ymean_by_xu{i},xu);            
    title(sprintf('%i %s %.3f %s',round(datadegmean(i)),'deg;',tunStrngth(i),'deg;'))
    hh(i) = vline(round(datadegmean(i)),':r');
    if i>1
        set(hh(1:end-1),'color','b'); set(hh(i),'color','r')
    end        
end

