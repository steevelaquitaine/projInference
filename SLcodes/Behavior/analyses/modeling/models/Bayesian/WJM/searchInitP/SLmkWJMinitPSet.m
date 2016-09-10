
%SLmkWJMinitPSet.m
%
% author: steeve laquitaine
%   date: 150715
%purpose: create reasonable combinations of kappa initial parameters to test
%         for WJM Bayes model
%         
%  usage:
%
%        initPset = SLmkWJMinitPSet(0:10:90,0:10:90,0:10:90)
%
%note: 
%        run with o = SLsearchWJMBayesinitP(initPset);

function initPset = SLmkWJMinitPSet(kappa24s,kappa12s,kappa06s)

%---------------------------------------------------------------
%Grid search of best initial kappa parameters for WJM Bayes model
%task conditions (ONCE)
%3 coh, 4 priors, motor noise, and random estimation
%[0.24 0.12 0.06 80 40 20 10 km prand]
%
%priors are assumed accurate (fixed)
%coherence is fixed
%prand is fixed at 0.0001
%only kappa and km change
%--------------------------------------------------------------
i = 0;
for kappa24 = kappa24s; %18 deg to 2 deg.
    for kappa12 = kappa12s; %18 deg to 2 deg.
        for kappa06 = kappa06s; %flat to 6 deg. 
            
            i = i + 1;      
            initPset(i,:) = [kappa24 kappa12 kappa06];
            
        end
    end
end
