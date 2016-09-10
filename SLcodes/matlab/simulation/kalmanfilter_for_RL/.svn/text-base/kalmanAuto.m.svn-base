function [negLogLike, parameters, exploration, td, exploitation] = kalmanAuto(reward, choice, numChoice, nodisplay)
%kalmanAuto  estimates parameters of Kalman filter model on the
%   basis of behavioral data from a group of participants. If data from a
%   group of participants is provided, the function assumes that all
%   parameters are the same across pariticipants except for explaration
%   parameter of Softmax. To fit separate parameters to each participant,
%   execute this function separately for each participant.
%
%   [NLL, P, Ep] = kalmanAuto(reward, choice, numChoice)
%   Input with behavioral data: reward, choice, numChoice
%       reward : a matrix of received rewards of participants
%       choice : a matrix of choices of participants
%           In the matrix, row represents trial and column is participant.
%           If a participant did not make a choice, the value of choice in
%           that trial becomes 0.
%       numChoice : the number of choices in the experiment
%       To see an example, please type "indirectActorAuto('example')".
%
%   returns negative log likelihood NLL, 
%   an array of parameters P:
%   Parameters : the array of parameters
%       parameters(1) is standard deviation of score,
%       parameters(2) is istandard deviation of diffusion,
%       parameters(3) is decaying rate of action values,
%       parameters(4) is converging value of action values,
%       parameters(5) is initial value of mean of score
%       and parameters(6) is initial value of standard deviation of score
%   ,and an array of exploration parameters Ep for all participants.
%
%   The function also generates information useful for the analysis of fMRI
%   data.
%
%   [NLL, P, Ep, T, Ex] = kalmanAuto(reward, choice, numChoice)
%   returns a matrix of TD(temporal difference) T
%   and a matrix of exploitations Ex.
%       In the matrix, row represents trial and column is participant.
%       exploitation is 1 when a participant exploits,
%       exploitation is 0 when a participant explores,
%       and exploitation is -1 when a participant does not make a choice
%
%   By default the function opens a figure with visualization of choice
%   probabilities predicted by the model and displays the parameters being
%   searched. To supress it, you can specify 4th optional parameter
%   'nodisplay'. indirectActorAuto(reward, choice, numChoice, nodisplay)
%   when nodisplay = 0, figure is generated (default)
%        nodisplay = 1, no figure is generated.
%
%   Jee Hoon, Yoo in University of Bristol, September 2008

if nargin == 1 & strcmp (reward, 'example')
 disp ('Let us consider the simple behavioral data. Only 2 participants take part');
 disp (' in the experiment, and in each trial they have to select one of 4 choices during ');
 disp (' the 60 trials. In this case, the matrix of reward is:');
 disp (' ');
 disp (' reward = [  72	130');
 disp ('            120	120');
 disp ('            122  64');
 disp ('            ...');
 disp ('            74  111 ]');
 disp (' ');
 disp (' Also, the matrix of choice is:');
 disp (' ');
 disp (' choice = [ 3 1');
 disp ('            4 2');
 disp ('            4 0');
 disp ('            ...');
 disp ('            1 2 ]');
 disp (' ');
 disp ('Note that the 2nd participant in 3rd trial did not make a choice. ');
 disp ('In this case, the script will ignore the reward in that trial');
 disp (' ');
 disp (' Lastly, the value of numChoice is:');
 disp (' ');
 disp (' numChoice = 4');
 disp (' ');
 disp (' There is an example data for test run. To see how it runs, type');
 disp (' in MATLAB: kalmanAuto(''testrun'');');

elseif nargin == 1 & strcmp (reward, 'testrun')
    load('example.mat');
    kalmanAuto(reward, choice, numChoice);
else
    
if nargin == 3
    nodisplay = 0; % set to default
end

disp ('Finding staring points of parameters...');
initValues = rand(10,6);
initResult = zeros(1, 10);
% use 10 different sets of parameters to find a starting point

maxReward = max(max(reward));
minReward = min(min(reward));
initRange = maxReward - minReward;
maxStd    = 10;
minStd    = 0;
StdRange  = maxStd - minStd;

for i = 1:10
    initValues(i,1) = initValues(i,1) * StdRange + minStd;
    initValues(i,2) = initValues(i,2) * StdRange + minStd;
    initValues(i,4) = initValues(i,4) * initRange + minReward;
    initValues(i,5) = initValues(i,5) * initRange + minReward;
    initValues(i,6) = initValues(i,6) * StdRange + minStd;
    initResult(i) = kalman(initValues(i,:), reward, choice, numChoice, nodisplay);
end
[minV minI] = min(initResult);

disp ('Initial parameters are:');
initParameters = initValues(minI,:)
% initial parameters are set to the ones with minimum negative log likelihood. 
    
disp ('Finding parameters in progress...');
parameters = fminsearch(@kalman, initParameters, [], reward, choice, numChoice, nodisplay);
% get the optimum parameters by using fminsearch.
disp ('Found parameters are:');
parameters

[negLogLike exploration] = kalman(parameters, reward, choice, numChoice, 1);
% get the best negative log likelihood and the values of exploration
% parameter by the obtained optimum parameters.
disp ('Found exploration parameters are:');
exploration

disp ('Found the best negative log likelihood is:');
negLogLike

[trials numOfData] = size(choice);
% get the number of trials and participants

td              = zeros(trials, numOfData);
exploitation    = zeros(trials, numOfData);

for i = 1:numOfData
    [indvResult, td(:, i), exploitation(:, i)] = kalmanIndv(exploration(i), reward(:, i), choice(:, i), numChoice, parameters);
    % get matricies of TD and exploitations by the optimum parameters.
end

if (nodisplay == 0)
    disp ('The represented graph is the example of the change of probability plot ');
    disp ('in 1st participant. To see how to use probplot, type:');
    disp ('resultplot(''example'');');

    resultplot(reward(:, 1), choice(:, 1), numChoice, parameters, exploration(1));
    % example of probplot
end

end
