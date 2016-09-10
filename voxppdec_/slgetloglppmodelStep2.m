

%slgetloglppmodelStep2.m
%
%
%
% author : steeve laquitaine
%purpose : get the logl of the pp model for parameters rho, tau and sigma


function [nglogl,fp] = slgetloglppmodelStep2(eta,nu,fp)

%get params. We substract 20 because we add 20 before entering the function
%to increase step size to 1 instead of actual parameter magnitude*5%. This
%make the search more explorative and speed up convergence
% fp = fp - 20;
% rho = fp(1);
% tau = fp(2:end-1);
% sigma = fp(end);
rho = fp(1) -  2*100;
tau = fp(2:end-1,1) - 1.4*100;
sigma = fp(end) - 1.8*100;

%sanity checks
if rho < 0 || rho >1 || any([rho;tau;sigma]<0)
    nglogl = inf;
    fprintf('%s \n', 'rho cannot be <0')
    return
end

%# of channels
K = size(eta,1);

%# of voxels
Nv = size(nu,1);

%# of instances
Ni = size(nu,2);

%probability of p_eta given the model parameters
p_eta = mvnpdf(eta',zeros(1,K),sigma^2*eye(K,K));

%probability of nu given the model parameters
SIGMA = rho*(tau*tau')+ (1-rho)*times(eye(Nv,Nv),tau*tau');

%get probability of each nu. The covariance matrix SIGMA must be 
%symmetric, positive definite so make sure SIGMA is a valid covariance
%matrix.
[~,err] = cholcov(SIGMA,0);
if err~=0
    fprintf('%s \n', 'Sigma was not a valid covariance matrix')
    nglogl = inf; return
end
p_nu = mvnpdf(nu',zeros(Ni,Nv),SIGMA);

%get the logl of the residual noises observed given the model params
%add different lapse rates lpss to deal with cases when data are never 
%predicted by the model (nglogl = - inf) due to e.g., numerical 
%precision. p_eta and p_nu do not have the same lapses because p_nu 
%is typically much lower (most values are near 0) that p_eta because 
%it is produced by a multivariate Gaussian on on many more dimensions 
%than p_eta so we re-scale to improve the fit of rho
% lps_eta = 0;
% lps_nu = 0;
% p_eta = p_eta + lps_eta;
% p_nu = p_nu + lps_nu;
% nglogl = - sum(log(p_eta)+log(p_nu)) + lps_eta;
nglogl = - sum(log([p_eta(:); p_nu(:)]));
% hold on; plot(nglogl,'o')
% drawnow

%objective value Blow-up
if nglogl==-inf
    dbstack
    keyboard
end

%print
fprintf('%s       %s        %s     %s     \n','nglogl','rho','sigma','mean(tau)')
fprintf('%.005f  %.006f   %.006f %.006f \n',nglogl,rho,sigma,nanmean(tau))