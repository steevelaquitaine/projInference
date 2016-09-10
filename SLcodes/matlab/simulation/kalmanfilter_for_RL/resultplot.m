function [] = resultplot(reward, choice, numChoice, parameters, exploration)
% resultplot plots the changes of probabilities and action values
%
% resultplot(reward, choice, numChoice, learningRate, exploration)
%   Inputs:
%       Behavioural data: reward, choice, numChoice
%           reward : the array of received rewards of a participant
%           choice : the array of choices of a participant
%           numChoice : the number of choices in the experiment
%
%       Model Parameters : parameters, exploration
%           parameters : other parameters from kalman function
%           exploration : obtained exploration parameter of Softmax method
%           from kalmanAuto
%
%   Outputs:
%       Two plots in each row represent the changes in probabilities predicted 
%       by the model and action values of one choice. For example, in the 1st row, 
%       there are two plots for 1st choice. 
%       In the graph of probability, the black line describes the changes
%       in probabilities, and the vertical line indicates that the participant 
%       chose the corresponding alternative in a given trial. In the graph of 
%       action value, the black line depicts the changes in action values, and 
%       the star represents the reward in each trial. For instance, the 
%       participant visualized above chose the 2nd alternative on 10th trial, 
%       hence two graphs corresponding to the 2 alternative (in 2nd row) show 
%       the vertical line and the star in 10th trial respectively.
%
%   Jee Hoon, Yoo in University of Bristol, September 2008

if nargin == 1 & strcmp (reward, 'example')
 disp ('Suppose that you have data for all participant; you probably have the matricies of');
 disp (' reward and choice, also the array of exploration parameters. If you want to');
 disp (' plot the change of probabilities and the change of action values in 5th participant,');
 disp (' type: resultplot(reward(:, 5), choice(:, 5), numChoice, learningRate, exploration(5));');

else
[negLogLike, td, exploitation, probs, m] = kalmanIndv(exploration, reward, choice, numChoice, parameters);
% get information from indirectActorIndv function

[trials numChoice] = size(probs);
% get the number of trials and choices.

row = numChoice;
col = 2;
% get the number of rows and columns.

for i = 1:numChoice
    probi = i*2 - 1;
    subplot(row, col, probi);
    plot(probs(:,i), 'k');
    hold on;
    
    for j = 1:trials
        if (choice(j) == i)
            plot(j, [0:0.01:1], 'b-.');
        end
    end
    hold off;
      
    if (i == 1)
        title('Probability');
    end      
    legend('probability', 'choice in trial');
    xlabel('trial');
    ylabel(['choice ' int2str(i)]);
    
    mi = i*2;
    subplot(row, col, mi);
    plot(m(:,i), 'k');
    hold on;
    
    for j = 1:trials
        if (choice(j) == i)
            plot(j, reward(j), '*');
        end
    end
    hold off;
      
    if (i == 1)
        title('Action value');
    end      
    legend('action value', 'reward in trial');
    xlabel('trial');
end

end