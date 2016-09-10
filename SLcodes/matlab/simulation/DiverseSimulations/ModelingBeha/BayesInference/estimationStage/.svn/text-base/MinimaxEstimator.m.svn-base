
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
    
function MinimaxEstimator(X)

% pdf.
s=tabulate(X);
p=s(:,3)/100;

% minimax estimator
x.space=unique(X);
for i=1:numel(x.space)
    
    % Loss function. Absolute error loss.
    loss(:,i)=abs(x.space-x.space(i));
    
    % Risk. Risk is defined as the average loss over x given parameter x(i)
    Risk(:,i)=sum(loss(:,i).*p);
end

% minimize Risk (over estimators)
Estimator=x.space(Risk==min(Risk));

% for pair dataset
Estimator=mean(Estimator);


% Talk
fprintf('\n %13s %12.8f', 'Median:',median(X))
fprintf('\n %13s %12.8f \n', 'Estimator:',Estimator)
    





