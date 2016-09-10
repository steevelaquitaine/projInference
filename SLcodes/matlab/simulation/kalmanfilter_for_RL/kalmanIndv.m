function [negLogLike, td, exploitation, probs, mRec] = kalmanIndv(exploration, reward, choice, numChoice, parameters)
%kalmanIndv  Negative log likelihood of Kalman filter model on the
%   basis of behavioral data from a participants. 
%
%   NLL = kalmanIndv(exploration, reward, choice, numChoice, parameters)
%   returns negative log likelihood NLL given behavioral data and parameters.
%
%   Behavioral data: reward, choice, numChoice
%       reward : the array of received rewards of a participant
%       choice : the array of choices of a participant
%       numChoice : the number of choices in the experiment
%
%   Parameters : exploration, 6 other parameters
%       exploration : exploration parameter of Softmax method
%       parameters(1) : standard deviation of score
%       parameters(2) : standard deviation of diffusion
%       parameters(3) : decaying rate of action values
%       parameters(4) : converging value of action values
%       parameters(5) : initial value of mean of score
%       parameters(6) : initial value of standard deviation of score
%       
%   [NLL, T] = kalmanIndv(...)
%   returns the array of TD(temporal difference) T.
%
%   [NLL, T, Ex] = kalmanIndv(...)
%   returns the array of exploitations Ex.
%       exploitation is 1 when a participant exploits,
%       exploitation is 0 when a participant explores,
%       and exploitation is -1 when a participant does not make a choice
%
%   [NLL, T, Ex, Ps] = kalmanIndv(...)
%   returns the matrix of probabilities of alternatives.
%       In the matrix, row represents trial and column is the number of choices.
%
%   Jee Hoon, Yoo in University of Bristol, September 2008

stdo        = parameters(1);
stdd        = parameters(2);
lambda      = parameters(3);
theta       = parameters(4);
initMean    = parameters(5);
initStd     = parameters(6);

trials      = length(choice);
% the number of trials

mean        = zeros(1, numChoice);
std         = zeros(1, numChoice);
% the mean and standard deviation of expected reward

probs       = zeros(trials, numChoice);
% probability of each choice in trials
mRec    = zeros(trials, numChoice);
% action values for recoding

td              = zeros(trials, 1);
% temporal difference

exploitation    = zeros(trials, 1);
% whether a participant exploits

%%%%% update procedure %%%%%
mean(:)     = lambda * initMean + (1 - lambda) * theta;
mRec(1, :)  = lambda * initMean + (1 - lambda) * theta;
std(:)      = lambda^2 * initStd + stdd^2;
% the mean and standard deviation of score set to initial value.

probs(1, :) = 1/numChoice;
% probabilities in first trial is initialized to have equal values.

if (exploration < 0 || stdo < 0 || stdd < 0 || lambda < 0 || initStd < 0)
    negLogLike = 10000;
    if (exploration < 0)
        negLogLike = negLogLike - exploration;
    end
    if (stdo < 0)
        negLogLike = negLogLike - stdo;
    end
    if (stdd < 0)
        negLogLike = negLogLike - stdd;
    end
    if (lambda < 0)
        negLogLike = negLogLike - lambda;
    end
    if (initStd < 0)
        negLogLike = negLogLike - initStd;
    end
else
for i = 1:(trials-1)
    if (choice(i) == 0)
        exploitation(i) = -1;
        % exploitation is set to -1 when a participant does not decide.
    else
        [maxValue maxIndex] = max(mean);
        if (choice(i) == maxIndex)
            exploitation(i) = 1;
        end
        % only when a participant chooses the score box of the highest
        % action value, exploitation is 1.
        
        kt              = std(choice(i)) / (std(choice(i)) + stdo^2);
        td(i)           = reward(i) - mean(choice(i));
        mean(choice(i)) = mean(choice(i)) + kt * td(i);
        std(choice(i))  = (1 - kt) * std(choice(i));
        % only chosen box's mean and standard deviation is updated.
    end
    
    probs(i+1,:) = softmax(mean, exploration);
    % probabilities are evaluated by Softmax method 
    % with exploration parameter.
    
    mean        = lambda * mean + (1 - lambda) * theta;
    std         = lambda^2 * std + stdd^2;
    mRec(i+1,:) = mean;
    % decay of the mean and the standard deviation
end

if (choice(trials) == 0)
    exploitation(trials) = -1;
else
    [maxValue maxIndex] = max(mean);
    if (choice(trials) == maxIndex)
        exploitation(trials) = 1;
    end
    
    kt                              = std(choice(trials)) / (std(choice(trials)) + stdo^2);
    td(trials)                      = reward(trials) - mean(choice(trials));
    mRec(trials, choice(trials))    = mRec(trials, choice(trials)) + kt * td(trials);
end
% the last updating procedure for exploitation and td.

%%%%% evaluation procedure %%%%%
negLogLike = 0;
for i = 1:trials
    if (choice(i) ~= 0)
        negLogLike = negLogLike - log(probs(i, choice(i)));
    end    
    % only when a participant makes a choice,
    % negative log likelihood is evaluated.
end

end