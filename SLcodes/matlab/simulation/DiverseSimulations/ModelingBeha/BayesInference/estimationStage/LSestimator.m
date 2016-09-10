
     % Author: Steeve Laquitaine
       % Date: 130710
    % Purpose: Calculate the minimax estimator of a dataset. The minimax
               %estimator is the median of the dataset.
  % reference: Figueiredo., Lecture notes on Bayesian Estimation and 
               %Classification, October 2004.
               %http://www.lx.it.pt/~mtf/learning/Bayes_lecture_notes.pdf

    % X is dataset
    
    % e.g.,
    % x=[1 2 4 4 4 4];
    % X=normrnd(25,5,50,1) 

    
function LSestimator(X)

% If X are not integers make them integers. If we don't do that data very
% rarely repeat and pdf is not gaussian but uniform.
if rem(X,1) ~= 0 ; % then its not integer
    fprintf('\n %s',['Careful: when data are not integers pdf may be]',...
        'flat due to sampling noise even if data were generated from a',...
        'Gaussian. This should not affect the result, but the flat pdf displayed'...
        '[may be confusing'])
end
x.space=unique(X);
s=tabulate(X);
[dummy,ints]=intersect(X,x.space);


% pdf
p=s(ints,3)/100;


% Least square estimation
% Try some estimators
for i=1:numel(x.space)
    
    % Loss function. Squared error loss.
    loss(:,i)=(x.space-x.space(i)).^2;
    
    % Risk. Risk is defined as the average loss over x given parameter x(i)
    Risk(:,i)=sum(loss(:,i).*p);
end

% minimize Risk (over estimators)
Estimator=x.space(Risk==min(Risk));

% for pair dataset
Estimator=mean(Estimator);


% Talk
fprintf('\n %13s %12.8f', 'Mean:',mean(X))
fprintf('\n %13s %12.8f \n', 'Estimator:',Estimator)


% Show
close all;
hold all
plot(x.space,p,'linewidth',2)
    