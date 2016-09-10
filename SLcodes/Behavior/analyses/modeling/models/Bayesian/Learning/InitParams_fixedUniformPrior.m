
%InitParams_fixedUniformPrior.m


%(1) Input parameters (best matching parameters)
k0(1,:) = [kllh    klearnt kcardinal fractRand motorNoise weightTail];
%(4) 8x stronger likelihoods & best matching priors
k0(2,:) = [kllh.*8 klearnt kcardinal fractRand motorNoise weightTail];
%(7) 8x weaker likelihood & true priors
k0(3,:) = [kllh./8 klearnt kcardinal fractRand motorNoise weightTail];
