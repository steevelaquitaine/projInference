%slWait.m
%
% author: steeve laquitaine
%purpose: wait n seconds before next action
%
%  usage:
%
%       slWait(2)

function slWait(nsec)

%hold n sec
t0 = tic;
t = 0;
while t < nsec
    t = toc(t0);
end