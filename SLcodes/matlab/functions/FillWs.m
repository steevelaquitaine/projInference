function [ output_args ] = FillWs( data )
%FILLWS Summary of this function goes here
%   Detailed explanation goes here

% opposite function to GrabWs: free fields contained
% in structure file and put them on workspace


name=fieldnames(data);
for i=1: numel(name)
    assignin('base',char(name(i)),data.(name{i}))
end

end
