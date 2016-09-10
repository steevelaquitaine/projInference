
%SLbinnedCumSum.m
%
% author: steeve laquitaine
%   date: 150721
%purpose: bin data and sum data within bin
%
% usage:
%
%          xbinned = SLbinnedCumSum(1:100,20,'index')
%          xbinned = SLbinnedCumSum(1:360,1:36:360,'index')
%          xbinned = SLbinnedCumSum(0:0.1:1,[zeros(1,9);0.1:1/9:1],'values')

function [xbinned,binStart, binEnd] = SLbinnedCumSum(x,mybin,type)

nx = length(x);

if strcmp(type,'index')
    %case bin
    if length(mybin) == 1
        
        tmp = 1 : mybin-1 : nx;
        binStart = tmp(1:end-1); %bin starts
        binEnd = tmp(2:end); %bin ends       
        nbins = length(binStart);

        for i = 1 : nbins
            
            xbinned(i) = sum(x(binStart(i) : binEnd(i)));
            
        end
        
        %case bin values
    elseif length(mybin) > 1 && size(mybin,1) == 1
        
        binStart = mybin(1:end-1); %bin starts
        binEnd = mybin(2:end); %bin ends
        nbins = length(binStart);
        
        for i = 1 : nbins
            
            if max(mybin) >  nx
                fprintf('(SLbinnedCumSum) Please make sure that max bin index <= length(x) \n')
                keyboard
            end
            
            xbinned(i) = sum(x(binStart(i) : binEnd(i)));

        end
        
        
        %case bin index start and end
    elseif length(mybin) > 1 && size(mybin,1) == 2
        binStart = mybin(1,:); %bin starts
        binEnd = mybin(2,:);   %bin ends
        nbins = size(mybin,2);
        
        for i = 1 : nbins
            
            xbinned(i) = sum(x(binStart(i) : binEnd(i)));
            
        end
    end
    
    
    
    %case bin value start and end values
elseif strcmp(type,'values')
    
    if length(mybin) > 1 && size(mybin,1) == 2
        
        binStart = mybin(1,:); %bin starts
        binEnd = mybin(2,:);   %bin ends
        nbins = size(mybin,2);
        
        for i = 1 : nbins
            
            myValsPos = SLfindBetween(x,binStart(i),binEnd(i),'includeBinStart');
            xbinned(i) = sum(x(myValsPos));
            
        end
    end
    
    
    
end








