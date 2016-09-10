 

     % Author: Steeve Laquitaine
       % Date: 2013/07/07
    % Purpose: Bayesian estimation with Minimax. 
% Description: Minimax decision rule aims to minimize minimum maximum losses.

function [estimator]=MinimaxReadOut(x) 

%% data
% x=[1 2 2];
x = [1 2 2]; 



% llh
stat=tabulate(x); p=stat(:,3)/100;



% plot(data)
close all
figure('color','w')
bar(unique(x),p,'edgecolor','w')
box off

% -----------------
% Method 1
% -----------------
% true mean (equivalent to sum([x.*p(x)))
Mean=mean(x);
% draw
hold on
plot([Mean Mean],[0 1],'-b','linewidth',2,'markersize',10)
text(Mean, -0.075, num2str(Mean))
publishAxis


% -----------------
% Method 2: Least square estimator
% -----------------
% give more resolution to x because calculating Risk is actually 
% calculating integrals over continuous loss function.
% xcont=x;%min(x):0.01:max(x);
xcont=interp1(x,1:0.01:max(x));
stat=tabulate(xcont);
p=stat(:,3)/100;

% Least square estimator. 
for i=1:numel(xcont)
    
    % loop over possible estimator
    es = xcont(i);
    
    % calulate loss function
    Loss(i,:) = (xcont - es).^2;
    
    % calculate Risk.Risk is averaged 
    % over all population values of x.
    Risk(:,i) = sum(p.*Loss(i,:)');
end

%Apply MLH decision rule to get best estimate.
[dummy, minmax]=min(Risk);
estimator=xcont(minmax);


% Draw
plot([estimator estimator],[0 1],'-r','linewidth',2,'markersize',10)
text(estimator, -0.1, ['es:' num2str(estimator)])


