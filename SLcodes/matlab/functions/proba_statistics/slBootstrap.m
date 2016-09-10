

%slBootstrap.m
%
%
% author: steeve laquitaine
%   date: 160427
%purpose: create nBoot new datasets by boostrapping {x,y} data pairs
%
%  usage:
%
%  e.g., 
%         x = 1:1:360; y = x+randn(length(x),1); 
%         [stat,x,y] = slBootstrap('stat=sllinefit(x,y,1:1:360,[1 0 0])',x,y,360,1)
%
%Inputs:
%       x: predictor/independent variable
%       y: predicted/dependent variable
%
%Outputs:
%       stats: product of an operation on x and y
%
%Description:
%       "myfun" requires three outputs : stat, x and y and an arbitrary
%       number of input including x and y


function [stat,x,y] = slBootstrap(myfun,x,y,ndata,nBoot)

%number of observations
nRows = length(x);
x_raw = SLmakeColumn(x);
y_raw = SLmakeColumn(y);

%sample {x,y} pairs at each bootstrap
stats = [];
for i = 1 : nBoot            
    Ix_sampled = randi([1 nRows],ndata,1);           
    x = x_raw(Ix_sampled);   
    y = y_raw(Ix_sampled);
    
    %evaluate function and output stats over {x,y}
    eval(myfun);    
    
    %backup
    stats{i} = stat;
    x_boot(:,i) = x;
    y_boot(:,i) = y;
end
stat = stats;
