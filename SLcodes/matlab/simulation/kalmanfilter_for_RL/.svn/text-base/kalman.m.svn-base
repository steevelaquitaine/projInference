function [negLogLike, exploration] = kalman(parameters, reward, choice, numChoice, nodisplay)
%kalman  Negative log likelihood of Kalman filter model on the
%   basis of behavioral data from a group of participants. It assume that all
%   parameters are the same across pariticipants except for explaration
%   parameter of Softmax.
%
%   NLL = kalman(parameters, reward, choice, numChoice)
%   returns negative log likelihood NLL given behavioral data and parameters.
%
%   Behavioral data: reward, choice, numChoice
%       reward : the matrix of received rewards of participants
%       choice : the matrix of choices of participants
%           In the matrix, row represents trial and column is participant.
%       numChoice : the number of choices in the experiment
%
%   Parameters : the array of parameters
%       stdo : parameters(1), standard deviation of score
%       stdd : parameters(2), standard deviation of diffusion
%       decayParameter : parameters(3), decaying rate of action values
%       decayCenter : parameters(4), converging value of action values
%       initMean : parameters(5), initial value of mean of score
%       initStd : parameters(6), initial value of standard deviation of score
%
%   [NLL, E] = kalman(...)
%   returns the array of exploration parameters E.
%       each value in the array is the optimum value of each participant.
%
%   Jee Hoon, Yoo in University of Bristol, September 2008

if nargin == 4
    nodisplay = 0; % set to default
end

[trials numOfData] = size(choice);
% get the number of trials and participants

initValues = [0.2 0.5 0.8];
initResult = zeros(1, 3);
% use 3 different exploration parameters to find a starting point

negLogLike = 0;
for i = 1:numOfData
    for j = 1:3
        initResult(j) = kalmanIndv(initValues(j), reward, choice, numChoice, parameters);
    end
    [minV minI]                     = min(initResult);
    [exploration(i) indvResult]     = fminsearch(@kalmanIndv, initValues(minI), [], reward(:, i), choice(:, i), numChoice, parameters);
    % get the optimum value of exploration parameter 
    % and its negative log likelihood
    
    negLogLike                      = negLogLike + indvResult;
end

if (nodisplay == 0)
    disp(['Parameters = ' num2str(parameters) ', NLL = ' num2str(negLogLike)]);
end
