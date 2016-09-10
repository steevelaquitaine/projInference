function [ regressors  ] = Qlearning(action,outcome, magnitude)

%% Optimize alpha and beta
params_RL = fminsearch(@(params_RL)-1 * FindlogL(action,outcome,magnitude, params_RL(1),params_RL(2)),[0.01 0.01],optimset('MaxFunEvals',100000,'MaxIter',100000));

%% QLEARNING
A=action;
R=outcome;
M=magnitude;

alpha   =params_RL(1);
beta    =params_RL(2);
qA        =0.5;                                                                             
qB        =0.5;
Q          =0.5;
p           =0.5;
error   = [];

for t=1:length(A)
    if A(t)==1                                                                           % chose A
        error (t)        =R   (t) -  qA(t);                                                    % PE Computation  
        qA      (t+1)   =qA(t) + alpha * error(t);                              % Action Value updating (RL algorithm, Me)% Action Value updating (RL algorithm)
        qB      (t+1)   =qB (t)  - alpha * error(t); 
        Q        (t+1)   =qA (t+1);                                                                          % Prediction of reward probability for experienced choice
    else                                                                                       % chose B
        error(t)        =R   (t) -  qB(t);
        qB     (t+1)   =qB (t) + alpha * error(t);                              % Other Action Value updating (RL algorithm)
        qA     (t+1)   =qA(t) - alpha * error(t);                              % Action Value updating (RL algorithm, Me)% Action Value updating (RL algorithm)
        Q       (t+1)   =qB (t+1);
    end
   p(t)=exp(qA(t)*M(t)/beta)/(exp(qA(t)*M(t)/beta)+exp(qB(t)*(1-M(t))/beta));                              % Softmax rule integrating reward magnitude
end

qA(end)=[];
qB(end)=[];
Q(end)=[];
    
hold on; plot([qA' qB' Q' error']);
regressors=[A R qA' qB' Q' error' p'];
