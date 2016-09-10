
%sltrainPPmodelStep2.m
%
% author: steeve laquitaine
%purpose: train the voxel probabilistic population code's noise parameters
%         rho, tau and sigma
%
%usage:
%
%
%       [rho,tau,sigma,nglogl] = sltrainPPmodelStep2(b,svec,f_k_s,W)
%
%inputs :
%
%           b : Ni instance x Nv voxels matrix
%        svec :
%       f_k_s :
%           W :

function [rho,tau,sigma,nglogl] = sltrainPPmodelStep2(b,svec,f_k_s,W)
tic
%# of channels
K = size(W,2);

%# of instances
Ni = length(svec);

%we just play with equation 1 from Van Bergen et al 2015 Nature
%Neuroscience to express it as a system of linear equation to be solved by 
%conjugate gradient method.
%The predicted Ni x Nv bold response matrix b_bred usually match pretty
%well the actual simulated matrices b with known parameters
%b_pred = (W*(f_k_s(:,svec)+eta) + nu)';
eta = nan(K,Ni);
for i = 1 : Ni
    
    %calculate residual
    B = b(i,:)' - W*f_k_s(:,svec(i));
    
    %We first get back eta: Ni x K channels matrix 
    %W must be a square matrix so we multiply 
    %everything by W' to fit with the conjugate gradient
    %method, Hopefully W'W is not ill-conditioned
    %W'*B = W'*W*eta + W'*nu;
    %check that W'W is not ill-conditioned
    if cond(W'*W) < 1/((max(size(W'*W))*eps)) == 0
        fprintf('%s \n','(sltrainPPmodelStep2) Coefficient matrix is ill-conditioned for conjugate gradient method')
    end
    eta(:,i) = conjgrad(W'*W,W'*B,1e-10);   
    
end
%We now get back nu
nu = b' - (W*f_k_s(:,svec) + W*eta);

%get pp model rho,tau and sigma best parameters
[rho,tau,sigma,nglogl] = slfitRhoTauSigma(eta,nu);
toc