
% author: steeve laquitaine
%   date: 13052015
%purpose: indicate when "esc" key is pressed.
%
%  usage:
%
%           k = SLgetKeyPress

function k = SLgetKeyPress

k=0;

while ~k
    k = waitforbuttonpress;
    if ~ double(get(gcf,'currentcharacter'))==27
        k = 0;
    end
end