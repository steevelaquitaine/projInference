function [ y ] = Bin( mat,bin )
%UNTITLED2 Summary of this function goes here
%  Detailed explanation goes here

% mat is the matrix to bin 
% bin is the bin size

    size_mat=size(mat,1);
    
    % bin size
    %bin=10;
    io=nan(size_mat/bin ,1);
    
    % create the bin
    for i=1:size_mat/bin      
        io(i,1)=1+i*bin;
        iend(i,1)=i*bin;           
    end
    io(end)=[];
    io(2:(size_mat/bin),1)=io;
    io(1)=1;   
    trl_wind=[io iend];  

    for i=1:size_mat/bin  
        y1=mat(trl_wind(i,1):trl_wind(i,2),:);
        y(1:size(y1(:),1),i)=y1(:);    
    end
  