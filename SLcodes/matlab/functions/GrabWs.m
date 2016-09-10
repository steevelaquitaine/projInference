function [ s ] = GrabWs( input_args )
% GRABWS Summary of this function goes here
%   Detailed explanation goes here
    % Collect all data from the workspace into 
    % a structure file
    % reference: http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matlab-workspace-onto-a-stack

% Grab workspace
w = evalin('caller', 'whos');
names = {w.name};
s = struct;
for i = 1:numel(w)
    s.(names{i}) = evalin('caller', names{i});
end

end

