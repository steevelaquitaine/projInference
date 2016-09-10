%SLgetLogL_WJMrotatedDirTest6.m
%
%
% author: steeve laquitaine
%   date: 150726 updated 150802
%purpose: calculate logl of given motion direction estimation
%         data given WJM model with the model estimates distribution predictions
%
%  usage: 
%
%       [loglAll,fitP,Pest] = SLgetLogL_WJMrotatedDirTest6(data,s,c,kappa,kappa_s,prand,km,fp,posC,Cond)
%
%
% outputs:
%
%           Pest : Ncond by 360 predicted estimates (1:1:360)            
%
% usually Ntrials = 5000
% not enough memory when 328 dots and 1e4 trials.

%get logl of data
function [loglAll,fitP,Pest] = SLgetLogL_WJMrotatedDirTest6(data,s,c,kappa,kappa_s,prand,km,fp,posC,Cond)

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

%first stim direction x (with trial-noise) - (with optimized circ_vmrnd)
x1_24     = slCirc_vmrnd2(single(su24(1)),single(fp(1)),Ntrials,Ndots);
x1_12     = slCirc_vmrnd2(single(su12(1)),single(fp(2)),Ntrials,Ndots);
x1_6      = slCirc_vmrnd2(single(su6(1)),single(fp(3)),Ntrials,Ndots);
x1        = [{x1_24} {x1_12} {x1_6}];
motn      = (2/pi/besseli(0,km,1)*exp(km*cos(motevec)-km))'; %motor noise
motn      = motn/sum(motn);
prior80   = double(bsxfun(@plus,M,2/pi/besseli(0,fp(4),1)*exp(fp(4) * cos(svec - 3.93) - fp(4))));%Mean of prior is 225 deg (3.93 rads),need high precision
prior40   = double(bsxfun(@plus,M,2/pi/besseli(0,fp(5),1)*exp(fp(5) * cos(svec - 3.93) - fp(5))));
prior20   = double(bsxfun(@plus,M,2/pi/besseli(0,fp(6),1)*exp(fp(6) * cos(svec - 3.93) - fp(6))));
prior10   = double(bsxfun(@plus,M,2/pi/besseli(0,fp(7),1)*exp(fp(7) * cos(svec - 3.93) - fp(7))));
Pest      = nan(360,4*sum([Nsu{:}]));

for i = 1 : length(coh)
   tail   = tailtmp(i);
   amp    = amp1tmp(i)/besseli(0,fp(i),1);
   Pesttmp= getPestEachCoh(Ntrials, Ndots, Ns, svec, su{i}, Nsu{i} ,x1{i} ,rotrad{i}, coh(i),fp(i),prior80,prior40,prior20,prior10,prand,motn);%Pest coh 24%
   Pest(:,w{i}) = Pesttmp;    %estimates densities (Ncond by Nest)
end

Pest      = bsxfun(@rdivide,Pest,sum(Pest));               %proba
Pest      = (1 - prand)*Pest + prand*(1/360);              %randomness
Pest      = SLcircConv(Pest,motn(:,ones(1,size(Pest,2),'int8')));%motor noise
Pest      = Pest(:,ix)';       %keep useful est densities

%sanity check
if size(Mcond,1)~=size(Cond,1) || ~isempty(setdiff(single(Mcond),single(Cond),'rows'))
    fprintf('Something"s wrong. The conditions associated with the model and the data do not match')
    keyboard
end

%get logl
logl = nan(nCon,1);
for Ci = 1 : nCon

    Pest(Ci,Pest(Ci,:) <= 0) = 10^-320;         %get rid of <0 due to convolution
    Pest(Ci,:) = Pest(Ci,:)./sum(Pest(Ci,:),2); %proba
    logl(Ci) = sum(log(Pest(Ci,data(posC==Ci))));

end
loglAll = sum(logl);
fitP    = fp;

%get Probability of estimate for an individual coherence
function Pest = getPestEachCoh(Ntrials,Ndots,Ns,svec, su, Nsu,x1,rotrad,c,kappa,prior80,prior40,prior20,prior10,prand,motn)

global tail
global amp
global predEsttmp
global pos
taili    = tail;
ampi     = amp;
predEst  = predEsttmp;
Pest_80  = nan(360,Nsu); %estimate density
Pest_40  = nan(360,Nsu);
Pest_20  = nan(360,Nsu);
Pest_10  = nan(360,Nsu);
coherent = (rand(Ntrials,Ndots) < c);                       % which dots are moving coherently
rndDots  = rand(Ntrials,Ndots)*2*pi;

%each unique stim dir this coh
parfor i = 1 : Nsu
    
    %x measurements
    x  = slRotateMatrixValues(x1,rotrad(i),'dispOff');
    x  = coherent .* x + (1-coherent) .* rndDots;  % Ntrials by Ndots
    x  = permute(repmat(x, [1 1 Ns]),[3 1 2]);     % to make dimensions match: Ns by Ntrials by Ndots
    %like per dot
    %This code produces von Mises densities vm(u,k) based on the
    %equation vm=exp(k.*cos(x-u))./(2*pi.*besseli(0,k)); The code works for any
    %value of k (but not for inf).The equation is adjusted because of the
    %following numerical issues: when k>700, vm is NaN because besseli(0,k) and
    %exp(k.*cos(x-u)) reach numerical limits. exp(k.*cos(x-u)-k) scales vm
    %without changing its shape. besseli(0,k,1)) does same. The adjusted
    %equation and the exact equation yield exact same results except that the
    %adjusted equation works for large k (>>700).
    like_per_dot = taili + ampi*exp(kappa * cos(bsxfun(@minus,double(x),svec))-kappa);
    like = prod(like_per_dot,3);

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
    postmeans_40 = atan2(ytemp_40,xtemp_40);  % 1 by Ntrials
    Pest_40(:,i) = hist(postmeans_40,predEst);%estimate distributions

    %prior 20
    %--------
    %0.01 s faster than bsxfun
    posterior_20 = like.*prior20;    
    xtemp_20     = sum(bsxfun(@times, cos(svec), posterior_20));
    ytemp_20     = sum(bsxfun(@times, sin(svec), posterior_20));
    postmeans_20 = atan2(ytemp_20,xtemp_20);  
    Pest_20(:,i) = hist(postmeans_20,predEst); %estimate distributions
    
    %prior 10
    %--------
    %0.01 s faster than bsxfun
    posterior_10 = like.*prior10;
    xtemp_10     = sum(bsxfun(@times, cos(svec), posterior_10));
    ytemp_10     = sum(bsxfun(@times, sin(svec), posterior_10));
    postmeans_10 = atan2(ytemp_10,xtemp_10); % 1 by Ntrials)
    Pest_10(:,i) = hist(postmeans_10,predEst);%estimate distributions
    
end

Pest  = [Pest_80(pos,:) Pest_40(pos,:) Pest_20(pos,:) Pest_10(pos,:)]; %Pestimate
