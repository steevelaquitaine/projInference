

%slsimPPchannels.m
%
%
%author : steeve laquitaine
%purpose: simulate 8 channels (half rectified cosines 
%         raised to power 5) over direction circular space 
%         1:1:360 degs with direction preferences
%         [45 90 135 180 225 270 315 360].
%
%usage :
%
%       [f_k_s, f_k_s_noisy] = slsimPPchannels(0.5)
%
%inputs:
%
%       sigma: standard deviation of the assumed 
%              shared Gaussian noise (deviation 
%              from mean response) between
%              channels with same direction preference.
%
%outputs :
%
%           pp.f_k_s : matrix of 8 channels tuning functions
%                      (K channels x S stimulus directions)
%                      
%     pp.f_k_s_noisy : f_k_s with shared Gaussian noise between
%                      channels with same direction preference.
%                      %note that tuning responses can be < 0 because
%                      of the noise
%           pp.means : means deviation (Gaussian noise) from channel mean response
%                      360 directions x K channels matrix of zeros
%       pp.covMatrix : diagonal square covariance matrix defining the 
%                      deviations from channel mean response where the
%                      diagonals are the variances sigma^2 of each
%                      channel response over trial instances
%                      8 channels x 8 channels matrix
%
%ref : equation 8, Ruben S van Bergen et al., 2015, Nature Neuroscience

function pp = slsimPPchannels(sigma)

%stimulus direction space
pp.s = 1:1:360;

%assumed 8 tuning functions
%with at least one at 180 deg
pp.phi_k = [45 90 135 180 225 270 315 360]';
pp.K = length(pp.phi_k);

%directions tuning functions 
%s - phi_k for each phi_k and s
phi_k = pp.phi_k;
s = pp.s;
sTophi_k = bsxfun(@(s,phi_k) s-phi_k,s,phi_k);
pp.f_k_s = max( 0,cos(pi*(sTophi_k)/180) ).^5;

%noise shared by channels with same tunings
%zero mean and all channels have same noise 
%std sigma^2;
%eta_k is a K channels x 360 directions
pp.means = zeros(360,pp.K);
pp.covMatrix = sigma^2*eye(pp.K,pp.K);
eta = mvnrnd(pp.means,pp.covMatrix)'; 

%noisy channel tuning functions
pp.f_k_s_noisy = pp.f_k_s + eta;
