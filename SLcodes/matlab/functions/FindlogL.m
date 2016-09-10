
%TEST1 Summary of this function goes here

%   Detailed explanation goes here

% load sub101209F
% load A_sub101209F
% load R_sub101209F

% e.g. [logL, regressors] = instr_RL(A_sub101209F,R_sub101209F,sub101209F(:,4)/100,alpha,beta);

function [log_L] = FindlogL(action,outcome,magnitude,alpha,beta)
%lOG_L, P, Qa,
r=outcome;
A=action;
M=magnitude;

log_L=0;
p=0.5;
qA=0.5;                                                                                                
qB=0.5;

% Alphas and betas, max log likelihood
for t=1:length(r);
    mag=M(t);
    if A(t)==1                                                               % correct choice, chose A
        error   =r (t)-qA;                                                    % PE Computation
        qA      =qA  + alpha * error;                                        % action Value updating (RL algorithm, Me)
        qB      =qB  - alpha * error;                                        % Shinsuke's Q learning
    else                                                                     % incorrect choice, chose B
        error   =r (t)-qB;
        qB      =qB  + alpha * error;                                        % alternative action Value updating (RL algorithm)
        qA      =qA  - alpha * error;                                        % Shinsuke's Q learning
    end
    Qa(t)=qA;
    Qb(t)=qB;
    p(t,1)=exp(qA*mag/beta)/(exp(qA*mag/beta)+exp(qB*(1-mag)/beta));              % Softmax rule integrating reward magnitude
end
log_L=log_L+mean(A.*log(p)+(1-A).*log(1-p));                             % update likelihood

end
