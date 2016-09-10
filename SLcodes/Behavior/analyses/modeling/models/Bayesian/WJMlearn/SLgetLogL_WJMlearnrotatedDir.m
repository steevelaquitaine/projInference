%SLgetLogL_WJMlearnrotatedDir.m
%
%
% author: steeve laquitaine
%   date: 150909
%purpose: calculate logl of given motion direction estimation
%         data given WJM model with the model estimates distribution predictions
%
%  usage:
%
%       [loglAll,fitP,Pest] = SLgetLogL_WJMlearnrotatedDir(data,s,c,kappa,kappa_s,prand,km,fp,posC,Cond)
%
%
% outputs:
%
%           Pest : Ncond by 360 predicted estimates (1:1:360)
%
% usually Ntrials = 5000
% not enough memory when 328 dots and 1e4 trials.

%get logl of data
function [neglogl,fitP,Pest] = SLgetLogL_WJMlearnrotatedDir(data,s,c,kappa,kappa_s,prand,km,fp,posC,Cond)

%Parameter check
%Get rid of iterations with unreasonable params.
%jump to next iteration.
%Constrain some params between 0 and 1, and some > 0.
%and case a parameter that is not cardinal prior strength is missing (NaN)
if prand > 1
    neglogl = 1e09;
    return
end
if any(fp < 0)
    neglogl = 1e09;
    return
end
if any(isnan(fp))
    neglogl = 1e09;
    fprintf('%s \n','(SLgetLogL_WJMlearnrotatedDir) A param is NaN')
    keyboard
end

global pmean
global su24
global su12
global su6
global Ntrials
global Ndots
global Ns
global nCon
global svec
global coh
global tail
global tailtmp
global amp
global amp1tmp
global rotrad
global motevec
global su
global Nsu
global w
global M
global ix
global Mcond
global shifttmp
global shift

%first stim direction x (with trial-noise) - (with optimized circ_vmrnd)
x1_24     = slCirc_vmrnd2(su24(1),fp(1),Ntrials,Ndots);
x1_12     = slCirc_vmrnd2(su12(1),fp(1),Ntrials,Ndots);
x1_6      = slCirc_vmrnd2(su6(1) ,fp(1),Ntrials,Ndots);
x1        = [{x1_24} {x1_12} {x1_6}];

%-------------------------------
%prior and motor representations
%-------------------------------
if km <= 1e300
    motn      = (2/pi/besseli(0,km,1)*exp(km*cos(motevec)-km))'; %motor noise
    motn      = motn/sum(motn);
else
    %delta density
    motn             = zeros(length(motevec),1);
    motn(motevec==0) = 1;
end

%prior 80
if fp(2) <= 1e300
    prior80   = bsxfun(@plus,M,2/pi/besseli(0,fp(2),1)*exp(fp(2) * cos(svec - pmean) - fp(2))); %Mean of prior is 225 deg (3.93 rads),need high precision
else
    %delta density
    prior80 = zeros(length(svec),1);
    prior80(svec==pmean) = 1;
    prior80 = prior80(:,ones(1,Ntrials));
end

%prior 40
if fp(3) <= 1e300
    prior40   = bsxfun(@plus,M,2/pi/besseli(0,fp(3),1)*exp(fp(3) * cos(svec - pmean) - fp(3)));
else
    %delta density
    prior40 = zeros(length(svec),1);
    prior40(svec==pmean) = 1;
    prior40 = prior40(:,ones(1,Ntrials));
end

%prior 20
if fp(4) <= 1e300
    prior20   = bsxfun(@plus,M,2/pi/besseli(0,fp(4),1)*exp(fp(4) * cos(svec - pmean) - fp(4)));
else
    %delta density
    prior20 = zeros(length(svec),1);
    prior20(svec==pmean) = 1;
    prior20 = prior20(:,ones(1,Ntrials));
end

%prior 10
if fp(5) <= 1e300
    prior10   = bsxfun(@plus,M,2/pi/besseli(0,fp(5),1)*exp(fp(5) * cos(svec - pmean) - fp(5)));
else
    %delta density
    prior10 = zeros(length(svec),1);
    prior10(svec==pmean) = 1;
    prior10 = prior10(:,ones(1,Ntrials));
end
Pest      = nan(360,4*sum([Nsu{:}]));

%each coherence
for i = 1 : length(coh)
    tail   = tailtmp(i);
    amp    = amp1tmp(i)/besseli(0,fp(1),1);
    shift  = shifttmp{i};
    Pesttmp = getPestEachCoh(Ntrials, Ndots, Ns, svec, su{i}, Nsu{i} ,x1{i} ,rotrad{i}, coh(i),fp(1),prior80,prior40,prior20,prior10,prand,motn);%Pest coh 24%
    Pest(:,w{i}) = Pesttmp;    %estimates densities (Ncond by Nest)
end

Pest      = bsxfun(@rdivide,Pest,sum(Pest));               %proba
Pest      = (1 - prand)*Pest + prand*(1/360);              %randomness
Pest      = SLcircConv(Pest,motn(:,ones(1,size(Pest,2)))); %motor noise
Pest      = Pest(:,ix)';                                   %keep useful est densities

%sanity check
if size(Mcond,1)~=size(Cond,1) || ~isempty(setdiff(Mcond,Cond,'rows'))
    fprintf('Something"s wrong. The conditions associated with the model and the data do not match')
    keyboard
end

%get logl
logltmp = nan(nCon,1);
for Ci = 1 : nCon
    
    Pest(Ci,Pest(Ci,:) <= 0) = 10^-320;         %get rid of <0 due to convolution
    Pest(Ci,:) = Pest(Ci,:)./sum(Pest(Ci,:),2); %proba
    logltmp(Ci) = sum(log(Pest(Ci,data(posC==Ci))));
    
end

neglogl = - sum(logltmp);
fitP = fp;

%get Probability of estimate for an individual coherence
function Pest = getPestEachCoh(Ntrials,Ndots,Ns,svec, su, Nsu,x1,rotrad,c,kappa,prior80,prior40,prior20,prior10,prand,motn)

global tail
global amp
global predEsttmp
global pos
global shift
taili    = tail;
ampi     = amp;
predEst  = predEsttmp;
shifti   = shift;
Pest_80  = nan(360,Nsu); %estimate density
Pest_40  = nan(360,Nsu);
Pest_20  = nan(360,Nsu);
Pest_10  = nan(360,Nsu);
coherent = (rand(Ntrials,Ndots) < c);  %which dots are moving coherently
rndDots  = rand(Ntrials,Ndots)*2*pi;

%x measurements
x  = x1;
x  = coherent .* x + (1-coherent) .* rndDots;  % Ntrials by Ndots
x  = permute(repmat(x, [1 1 Ns]),[3 1 2]);     % to make dimensions match: Ns by Ntrials by Ndots

keyboard

%This code produces von Mises densities vm(u,k) based on the
%equation vm=exp(k.*cos(x-u))./(2*pi.*besseli(0,k)); The code works for any
%value of k (but not for inf).The equation is adjusted because of the
%following numerical issues: when k>700, vm is NaN because besseli(0,k) and
%exp(k.*cos(x-u)) reach numerical limits. exp(k.*cos(x-u)-k) scales vm
%without changing its shape. besseli(0,k,1)) does same. The adjusted
%equation and the exact equation yield exact same results except that the
%adjusted equation works for large k (>>700).
like_per_dot = taili + ampi*exp(kappa * cos(bsxfun(@minus,x,svec))-kappa);
like1        = prod(like_per_dot,3);

%each displayed direction
parfor i = 1 : Nsu
    
    %Like Matx of other disp.dir. are circ shifts of the 1st dir like
    like = circshift(like1,[-shifti(i) 0]);
    
    %prior 80
    %--------
    %0.01 s faster than bsxfun
    posterior_80 = like.*prior80;
    xtemp_80     = sum(bsxfun(@times, cos(svec), posterior_80));
    ytemp_80     = sum(bsxfun(@times, sin(svec), posterior_80));
    postmeans_80 = atan2(ytemp_80,xtemp_80); % 1 by Ntrials)
    Pest_80(:,i) = hist(postmeans_80,predEst);     %estimate distributions
    
    %prior 40
    %--------
    %0.01 s faster than bsxfun
    posterior_40 = like.*prior40;
    xtemp_40     = sum(bsxfun(@times, cos(svec), posterior_40));
    ytemp_40     = sum(bsxfun(@times, sin(svec), posterior_40));
    postmeans_40 = atan2(ytemp_40,xtemp_40);  
    Pest_40(:,i) = hist(postmeans_40,predEst);
    
    %prior 20
    %--------
    %0.01 s faster than bsxfun
    posterior_20 = like.*prior20;
    xtemp_20     = sum(bsxfun(@times, cos(svec), posterior_20));
    ytemp_20     = sum(bsxfun(@times, sin(svec), posterior_20));
    postmeans_20 = atan2(ytemp_20,xtemp_20);
    Pest_20(:,i) = hist(postmeans_20,predEst);
    
    %prior 10
    %--------
    %0.01 s faster than bsxfun
    posterior_10 = like.*prior10;
    xtemp_10     = sum(bsxfun(@times, cos(svec), posterior_10));
    ytemp_10     = sum(bsxfun(@times, sin(svec), posterior_10));
    postmeans_10 = atan2(ytemp_10,xtemp_10); 
    Pest_10(:,i) = hist(postmeans_10,predEst);
    
end

Pest  = [Pest_80(pos,:) Pest_40(pos,:) Pest_20(pos,:) Pest_10(pos,:)]; %Pestimate
