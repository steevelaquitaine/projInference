

%norm increases as rho (shared voxels variability) increases 


%random cov
Omega0 = rand(100,100);
parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega0))))*exp(-0.5*(b(:,i) - mu_i)'*Omega0*(b(:,i) - mu_i));    
end
nglogl0 = - sum(p_i);


parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega))))*exp(-0.5*(b(:,i) - mu_i)'*Omega*(b(:,i) - mu_i));    
end
nglogl0 = - sum(p_i);


sigma = 0
Omega1 = rho*tau*tau' +  (1-rho)*times(eye(Nv,Nv),tau*tau') + sigma^2.*W*W';
parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega1))))*exp(-0.5*(b(:,i) - mu_i)'*Omega1*(b(:,i) - mu_i));    
end
nglogl1 = - sum(p_i);


sigma = 0.1
Omega2 = rho*tau*tau' +  (1-rho)*times(eye(Nv,Nv),tau*tau') + sigma^2.*W*W';
parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega2))))*exp(-0.5*(b(:,i) - mu_i)'*Omega2*(b(:,i) - mu_i));    
end
nglogl2 = - sum(p_i);

sigma = 0.2
Omega3 = rho*tau*tau' +  (1-rho)*times(eye(Nv,Nv),tau*tau') + sigma^2.*W*W';
parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega3))))*exp(-0.5*(b(:,i) - mu_i)'*Omega3*(b(:,i) - mu_i));    
end
nglogl3 = - sum(p_i);


sigma = 0.3
Omega4 = rho*tau*tau' +  (1-rho)*times(eye(Nv,Nv),tau*tau') + sigma^2.*W*W';
parfor i = 1 : length(svec)    
    mu_i = W*f_k_s(:,svec(i));        
    p_i(i) = (1/(sqrt(2*pi*det(Omega4))))*exp(-0.5*(b(:,i) - mu_i)'*Omega4*(b(:,i) - mu_i));    
end
nglogl4 = - sum(p_i);