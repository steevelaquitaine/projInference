
%SLgetLogL_WJM.m
%
%
% author: steeve laquitaine
%   date: 150726 updated 150802
%purpose: calculate logl of given motion direction estimation
%         data given WJM model with the model estimates distribution predictions
%
%  usage: 
%
%       [loglAll,fitP,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fp,posC)
%
%
% outputs:
%
%           PestGivModel : Ncond by 360 predicted estimates (1:1:360)            
%
% usually Ntrials = 5000
% not enough memory when 328 dots and 1e4 trials.

%get logl of data
function [loglAll,fitP,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fp,posC)

tic

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
parfor Ci = 1 : nCon
        
    %measurements
    coherent = (rand(Ntrials,Ndots) < c(Ci)); %which dots are moving coherently
    x = circ_vmrnd(s(Ci)*ones(Ntrials,Ndots),kappa(Ci));
    x = coherent .* x + (1-coherent) .* rand(Ntrials,Ndots)*2*pi; % Ntrials by Ndots
    
    %prior
    Mean_s = SLde2r(225,0);          %known by subj (fixed)
    prior = 2/pi/besseli(0,kappa(Ci),1)*exp(kappa_s(Ci) * cos(svec - Mean_s) - kappa_s(Ci));
    
    %Likelihoods L(s) = p(x|s)
    x = permute(repmat(x, [1 1 Ns]),[3 1 2]); % to make dimensions match: Ns by Ntrials by Ndots

    %like per dot
    %This code produces von Mises densities vm(u,k) based on the
    %equation vm=exp(k.*cos(x-u))./(2*pi.*besseli(0,k)); The code works for any
    %value of k (but not for inf).The equation is adjusted because of the
    %following numerical issues: when k>700, vm is NaN because besseli(0,k) and
    %exp(k.*cos(x-u)) reach numerical limits. exp(k.*cos(x-u)-k) scales vm
    %without changing its shape. besseli(0,k,1)) does same. The adjusted
    %equation and the exact equation yield exact same results except that the
    %adjusted equation works for large k (>>700).
    like_per_dot = (1-c(Ci))/2/pi + c(Ci)/2/pi/besseli(0,kappa(Ci),1) * exp(kappa(Ci) * cos(bsxfun(@minus,x,svec))-kappa(Ci)); %(1-c)/2pi + c*V(svec,x,kappa)
    like = prod(like_per_dot,3);
    
    %Posteriors p(s|x)
    posterior = bsxfun(@times, prior, like); % Ns by Ntrials
    
    %Estimation: posterior means (circular)
    xtemp = sum(bsxfun(@times, cos(svec), posterior));
    ytemp = sum(bsxfun(@times, sin(svec), posterior));
    posteriormeans = atan2(ytemp,xtemp); % 1 by Ntrials)
    posteriormeans = SLra2d(posteriormeans);
    
    %estimate distributions
    PestGivModel(Ci,:) = hist(posteriormeans,predEst);
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    PestGivModel(Ci,:) = (1-prand)*PestGivModel(Ci,:) + prand*(1/360);               %add random estimation
    PestGivModel(Ci,:) = SLcircConv(PestGivModel(Ci,:),vmPdfs(0:1:359,0,km,'norm')');%add motor noise

    fprintf('%s %i \n','Cond: ',Ci)
end

toc
keyboard

%get logl
for Ci = 1 : nCon
    PestGivModel(Ci,PestGivModel(Ci,:)<=0) = 10^-320;                   %get rid of negatives due to convolution
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    logl(Ci) = sum(log(PestGivModel(Ci,data(posC==Ci))));
end

loglAll = sum(logl);
fitP = fp; %best fitp



