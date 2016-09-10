

%slgetloglPPmodel
%
%
%usage:
%
%       p_i = slgetloglPPmodel(b,svec,f_k_s,W,[rho tau sigma]);


function [nglogl,p_i] = slgetloglPPmodel(b,svec,f_k_s,W,fp)

%fit parameters
rho = fp(1);
tau = fp(2:end-1);
sigma = fp(3);

%constrain fp > 0
if any(fp<0)
    fprintf('%s \n','Fit parameters should not be < 0. Omega must be >0 definite')
    return
end

%# of voxels
Nv = size(W,1);

%covariance matrix Omega
Omega = rho*tau*tau' +  (1-rho)*times(eye(Nv,Nv),tau*tau') + sigma^2*W*W';

%get likelihood of each trial instance
%bold response pattern given the generative model
b_disc = round(b*10)/10;
bold_space = min(b_disc (:)) : 0.1 : max(b_disc (:));
parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega))))*exp(-0.5*(b(:,i) - mu_i)'*Omega*(b(:,i) - mu_i));    
end

%calculate negative loglikelihood to minimize
nglogl = - sum(p_i);


%print prms and obj
fprintf('%i %.02f %.02f \n',nglogl,rho,sigma)
