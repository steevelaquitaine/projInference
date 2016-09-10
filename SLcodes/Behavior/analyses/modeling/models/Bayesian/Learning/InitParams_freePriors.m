

%purpose : set initial parameters for learning Bayesian model
%
% 12 parameters : 
% 3 kllh, 4 klearnt, 1 kcardinal, 1 fractRand, 1 motorNoise, 1 weightTail, 1 tau        


%Sets of initial parameters: best matching parameters are found by
%matching raw data with simulations graphically (called "Kreal").
%We pre-generate 27 reasonable sets of initial parameters.
%+-----+-------+----------------------+
%|     |       |          prior       |
%.-----|-------+-------+-------+------+
%|     |       | true  | strong| weak |
%.-----|-------+-------+-------+------+
%|     |true   |(1)t-t |(2)t-s |(3)t-w|
%.     |-------+-------+-------+------+
%| llh |strong |(4)s-t |(5)s-s |(6)s-w|
%|     |-------+-------+-------+------+
%|     |weak   |(7)w-t |(8)w-s |(9)w-w|
%+-----+-----------------------+------+

%set 10: likelihood and prior strengths are all same.
%eventually we may use later
%10 sets.

%note: strong priors (and weak priors) are 8 times stronger than best
%matching priors. 8x is the factor for which I see clear deviation
%of simulation from data.

%We used one intial value for probability of random estimation, motor
%noise and cardinal prior strength, that we think are relatively small
%values (high k for motor noise is low motor noise) by looking at the
%data.

%Input initial fit parameters (typically best matching initial values)
%check that all initial parameters are input
if length(initp)~=12
    fprintf('%s \n',['(slfitMaxLLH) Initial parameters are missing...',...
        'parameters for learning Bayesian model is:', '\n 3 kllh', '\n 4 klearnt','\n 1 kcardinal',...
        '\n 1 fractRand','\n 1 motorNoise','\n 1 weighTail','\n 1 learning tau'])
    keyboard
end
kllh       = initp(1:3);
klearnt    = initp(4:7);
kcardinal  = initp(8);
fractRand  = initp(9);
motorNoise = initp(10);
weightTail = initp(11);
tau        = initp(12);

%(1) Input parameters (best matching parameters)
k0(1,:) = [kllh    klearnt      kcardinal fractRand motorNoise weightTail tau];
%(2) 8x stronger priors & best matching llh
k0(2,:) = [kllh    8*klearnt    kcardinal fractRand motorNoise weightTail tau];
%(3) 8x weaker priors & best matching llh
k0(3,:) = [kllh    klearnt./8   kcardinal fractRand motorNoise weightTail tau];
%(4) 8x stronger likelihoods & best matching priors
k0(4,:) = [kllh.*8 klearnt      kcardinal fractRand motorNoise weightTail tau];
%(5) 8x stronger likelihoods & stronger priors
k0(5,:) = [[kllh   klearnt].*8  kcardinal fractRand motorNoise weightTail tau];
%(6) 8x stronger likelihoods & weaker priors
k0(6,:) = [kllh.*8 klearnt./8   kcardinal fractRand motorNoise weightTail tau];
%(7) 8x weaker likelihood & true priors
k0(7,:) = [kllh./8 klearnt      kcardinal fractRand motorNoise weightTail tau];
%(8) 8x weaker likelihood & stronger priors
k0(8,:) = [kllh./8 klearnt.*8   kcardinal fractRand motorNoise weightTail tau];
%(9) 8x weaker likelihood & weaker priors
k0(9,:) = [[kllh   klearnt]./8  kcardinal fractRand motorNoise weightTail tau];
%(10)likelihoods & priors are same.
k0(10,:) = [repmat(nanmean(k0(1,1:7)),1,7) kcardinal fractRand motorNoise weightTail tau];

output.initp = k0;
