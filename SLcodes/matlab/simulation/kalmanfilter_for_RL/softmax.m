function probs = softmax(actionValues, b)
%softmax  probability evaluation based on action values
%
%   P = softmax(actionValues, b)
%   returns the array of probabilities P given action values and
%       exploration parameter b.
%       the value of exploration parameter should be bigger than 0 and smaller
%       than 1.
%
%   Jee Hoon, Yoo in University of Bristol, September 2008

n = length(actionValues);
% the number of choices

total = 0;
for i = 1:n
    total = total + exp(b * actionValues(i));
    % get the sum of all action values.
end

for i = 1:n
    probs(i) = exp(b * actionValues(i)) / total;
    % evaluate each probability of action value.
end
