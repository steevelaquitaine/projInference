
%slPrintfStr.m
%
%
% author: steeve laquitaine
%purpose: print basic text
%
%usage:
%
%       slPrintfStr('slMean','Calculating mean')

function slPrintfStr(myfun,mytxt)

fprintf('%s \n',['(',myfun,')',mytxt])


