
%slsimPPCinstances.m
%
%
% author : steeve laquitaine
%purpose : generates S response instances x N voxels
%          for a vector of S stimulus directions
%          stimulus directions values must belong to 1:1:360 deg
%
%usage :
%
%         svec = [15*ones(10,1); 85*ones(10,1)];
%         W = rand(100,8); 
%         sigma = 0.5;
%         tau = rand(100,1); 
%         rho = 0.5;
%         [b_wres,b,pp] = slsimPPCinstances(svec,W,sigma,tau,rho)
%
%inputs
%   
%          svec : vector of trial instances of stimulus directions
%             W : a matrix of voxels x channel positive weights
%         sigma : standard deviation of the assumed 
%                 shared Gaussian noise (deviation 
%                 from mean response) between
%                 channels with same direction preference.
%
%
%outputs:
%
%             
%        b_wres : a matrix of S noisy response instances x N voxels 
%                 include noise eta shared between channels with same 
%                 direction preferences and residual noise that 
%                 is the deviation noises specific to individual voxels
%                 (with std tau_i for voxel i), as well as noise shared 
%                 globally among voxels irrespective of their tuning.
%             C : a matrix of the channel average responses to the
%                 vector of stimulus directions
%                 Ni instances x 8 channels
%             b : a matrix of S response instances x N voxels 
%            pp : structure containing model architecture and 
%                 parameters
%
%ref : Ruben S van Bergen et al., 2015, Nature Neuroscience

function [b_wres,C,b,pp] = slsimPPCinstances(svec,W,sigma,tau,rho)


%8 direction tuned channels
%note that tuning responses can be < 0 because
%of the noise
pp = slsimPPchannels(sigma);
f_k_s_noisy = pp.f_k_s_noisy';

%average channel responses to vector of stimulus directions
%Ni instances x 8 channels
C = pp.f_k_s(:,svec)';

%simulate a matrix b of Ni trial stimulus instances x N voxels 
Ni = length(svec);
N = size(W,1);
b = nan(Ni,N);
for i = 1 :  N   
    Wi = W(i,:);
    b(:,i) = sum(Wi(ones(1,Ni),:).*f_k_s_noisy(svec,:),2);
end

%Add S instances x N voxel noise matrix nu to the channel mean responses
%noise is sampled out of a multivariate Gaussian with mean 0 and 
%covariance SIGMA modelled by parameters rho and tau.
%The noise is composed of noise specific to individual voxels (with std
%tau_i for voxel i), and of noise shared among all voxels 
%irrespective of their tuning.
%
%Warning : need to be careful and use elementwise hadamard product 
%"times(eye(N,N),tau*tau')". The higher rho is and the more clustered by 
%stimulus directions the covariance is between pairs of voxels. When rho is
%low covariance is more scattered.
pp.meansres = zeros(Ni,N);
pp.covMatrixres = rho*(tau*tau') + (1-rho)*times(eye(N,N),tau*tau');
nu = mvnrnd(pp.meansres,pp.covMatrixres); 
b_wres = b + nu;

%save params
pp.rho = rho;
pp.sigma = sigma;
pp.tau = tau;
