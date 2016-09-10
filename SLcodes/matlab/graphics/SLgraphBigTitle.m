
%SLgraphBigTitle.m
%
%  author: steeve laquitaine
%    date: 2015/04/22

function SLgraphBigTitle(myText)

annotation(gcf,'textbox',...
    [0.45 0.95 0.15 0.05],...
    'String',{myText},'FitBoxToText','on','EdgeColor','none');