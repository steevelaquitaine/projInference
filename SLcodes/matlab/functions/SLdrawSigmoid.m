

% author: steeve laquitaine
%   date: 140815
%purpose: plot divisive or sigmoid function
%
%  usage: 
%
%       P1 = SLdrawSigmoid(50,linspace(1,100),'divisive')
%       P1 = SLdrawSigmoid(50,linspace(1,100),'sigmoid')

function  P2 = SLdrawSigmoid(k1,k2,varargin)

numk2 = numel(k2);

%divisive competition function
if sum(strcmp(varargin,'divisive')) == 1
    for i = 1 : numk2
        P2(i) = k2(i)./(k1+k2(i));
    end
    set(gcf,'color','w')
    plot(k2,P2,'r','linesmoothing','on')
    box offs
    axis square
    ylim([0 1])
end

%sigmoid activation function
if sum(strcmp(varargin,'sigmoid')) == 1
    for i = 1 : numk2
        P2(i) = exp(k2(i)) ./ (exp(k1) + exp(k2(i)));
    end
    set(gcf,'color','w')
    plot(k2,P2,'r','linesmoothing','on')
    box off
    axis square
    ylim([0 1])
end
