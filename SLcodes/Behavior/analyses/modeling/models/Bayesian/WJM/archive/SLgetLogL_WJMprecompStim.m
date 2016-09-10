
%SLgetLogL_WJMprecompStim.m
%
%
% author: steeve laquitaine
%   date: 150726 updated 150802
%purpose: calculate logl of given motion direction estimation
%         data given WJM model with the model estimates distribution predictions
%
%  usage: 
%
%       [loglAll,fitP,PestGivModel] = SLgetLogL_WJMprecompStim(data,s,c,kappa,kappa_s,prand,km,fp,posC)
%
%
% outputs:
%
%           PestGivModel : Ncond by 360 predicted estimates (1:1:360)            
%
% usually Ntrials = 5000
% not enough memory when 328 dots and 1e4 trials.

%get logl of data
function [loglAll,fitP,PestGivModel] = SLgetLogL_WJMprecompStim(data,s,c,kappa,kappa_s,prand,km,fp,posC)

%Simulating measurements x
Ndots = 328;      %nb of dots
Ntrials = 5000;   %nb of trials

%prior parameters
Ns = 100;
svec = linspace(1,360,Ns+1)';
svec = SLde2r(svec(1:end-1),0);  %vector of hypothesized motion directions

%predicted estimates space
predEst = 1:1:360;
nCon = max(posC);

%each condition
for Ci = 1 : nCon
        
    tic    
    %measurements
    coherent = (rand(Ntrials,Ndots) < c(Ci)); %which dots are moving coherently
    toci(1,Ci) = toc ;
    
    tic    
    x = circ_vmrnd(s(Ci)*ones(Ntrials,Ndots),kappa(Ci));
    toci(2,Ci) = toc ;
    
    tic    
    x = coherent .* x + (1-coherent) .* rand(Ntrials,Ndots)*2*pi; % Ntrials by Ndots
    toci(3,Ci) = toc ;
    
    %prior
    tic    
    Mean_s = SLde2r(225,0);          %known by subj (fixed)
    prior = 2/pi/besseli(0,kappa(Ci),1)*exp(kappa_s(Ci) * cos(svec - Mean_s) - kappa_s(Ci));
    toci(4,Ci) = toc;

    %Likelihoods L(s) = p(x|s)
    tic    
    x = permute(repmat(x, [1 1 Ns]),[3 1 2]); % to make dimensions match: Ns by Ntrials by Ndots
    toci(5,Ci) = toc ;

    %like per dot
    %This code produces von Mises densities vm(u,k) based on the
    %equation vm=exp(k.*cos(x-u))./(2*pi.*besseli(0,k)); The code works for any
    %value of k (but not for inf).The equation is adjusted because of the
    %following numerical issues: when k>700, vm is NaN because besseli(0,k) and
    %exp(k.*cos(x-u)) reach numerical limits. exp(k.*cos(x-u)-k) scales vm
    %without changing its shape. besseli(0,k,1)) does same. The adjusted
    %equation and the exact equation yield exact same results except that the
    %adjusted equation works for large k (>>700).
    tic
    like_per_dot = (1-c(Ci))/2/pi + c(Ci)/2/pi/besseli(0,kappa(Ci),1) * exp(kappa(Ci) * cos(bsxfun(@minus,x,svec))-kappa(Ci)); %(1-c)/2pi + c*V(svec,x,kappa)
    toci(6,Ci) = toc;

    tic
    like = prod(like_per_dot,3);
    toci(7,Ci) = toc ;
    
    %Posteriors p(s|x)
    tic
    posterior = bsxfun(@times, prior, like); % Ns by Ntrials
    toci(8,Ci) = toc; 

    %Estimation: posterior means (circular)
    tic
    xtemp = sum(bsxfun(@times, cos(svec), posterior));
    toci(9,Ci) = toc ;

    tic
    ytemp = sum(bsxfun(@times, sin(svec), posterior));
    toci(10,Ci) = toc ;

    tic
    posteriormeans = atan2(ytemp,xtemp); % 1 by Ntrials)
    toci(11,Ci) = toc ;

    tic
    posteriormeans = SLra2d(posteriormeans);
    toci(12,Ci) = toc ;

    %estimate distributions
    tic
    PestGivModel(Ci,:) = hist(posteriormeans,predEst);
    toci(13,Ci) = toc ;

    tic
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    toci(14,Ci) = toc ;

    tic
    PestGivModel(Ci,:) = (1-prand)*PestGivModel(Ci,:) + prand*(1/360);               %add random estimation
    toci(15,Ci) = toc ;

    tic
    PestGivModel(Ci,:) = SLcircConv(PestGivModel(Ci,:),vmPdfs(0:1:359,0,km,'norm')');%add motor noise
    toci(16,Ci) = toc ;

    tic
    fprintf('%s %i \n','Cond: ',Ci)
    toci(17,Ci) = toc ;
end

keyboard

%get logl
for Ci = 1 : nCon
    PestGivModel(Ci,PestGivModel(Ci,:)<=0) = 10^-320;                   %get rid of negatives due to convolution
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    logl(Ci) = sum(log(PestGivModel(Ci,data(posC==Ci))));
end

loglAll = sum(logl);
fitP = fp; %best fitp



