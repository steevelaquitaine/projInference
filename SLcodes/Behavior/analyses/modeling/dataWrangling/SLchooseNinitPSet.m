
%SLchooseNinitPSet.m
%
% author: steeve laquitaine
%   date: 150715
%purpose: choose N initial parameter sets among all combinations of kappa 
%         initial parameters to test for WJM Bayes model
%         
%  usage:
%
%        initPsetAll = SLmkWJMinitPSet(0:10:90,0:10:90,0:10:90)
%        initPset = SLchooseNinitPSet(initPsetAll,1:200)
%

function initPset = SLchooseNinitPSet(initPsetAll,Nset)

initPset = initPsetAll(Nset,:);