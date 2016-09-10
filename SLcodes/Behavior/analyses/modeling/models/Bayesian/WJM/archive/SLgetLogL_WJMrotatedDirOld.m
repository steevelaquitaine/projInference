%SLgetLogL_WJMrotatedDirOld.m
%
%
% author: steeve laquitaine
%   date: 150726 updated 150802
%purpose: calculate logl of given motion direction estimation
%         data given WJM model with the model estimates distribution predictions
%
%  usage: 
%
%       [loglAll,fitP,Pest] = SLgetLogL_WJMrotatedDirOld(data,s,c,kappa,kappa_s,prand,km,fp,posC)
%
%
% outputs:
%
%           Pest : Ncond by 360 predicted estimates (1:1:360)            
%
% usually Ntrials = 5000
% not enough memory when 328 dots and 1e4 trials.

%get logl of data
function [loglAll,fitP,Pest] = SLgetLogL_WJMrotatedDirOld(data,s,c,kappa,kappa_s,prand,km,fp,posC,Cond)


%Simulating measurements x
%Ndots = 328;      %nb of dots
%Ntrials = 5000;   %nb of trials

%Simulating measurements x
Ndots   = 328 ;  %nb of dots
Ntrials = 250;   %nb of trials

%prior parameters
Ns = 100;
svecdeg = linspace(1,360,Ns+1)'; 
svecdeg = svecdeg(1:end-1);
svec = SLde2r(svecdeg,0);  %vector of hypothesized motion directions

%predicted estimates space
nCon = max(posC);

%------
%priors
%------
Mean_s = SLde2r(225,0); % assumed known by subj
prior80 = 2/pi/besseli(0,fp(4),1)*exp(fp(4) * cos(svec - Mean_s) - fp(4));
prior40 = 2/pi/besseli(0,fp(5),1)*exp(fp(5) * cos(svec - Mean_s) - fp(5));
prior20 = 2/pi/besseli(0,fp(6),1)*exp(fp(6) * cos(svec - Mean_s) - fp(6));
prior10 = 2/pi/besseli(0,fp(7),1)*exp(fp(7) * cos(svec - Mean_s) - fp(7));
    

%--------------------------------------
%x for each coherence & directionss
%--------------------------------------
%calculate x (Ntrials by Ndots) for condition 0.24 & s(1)
%This preserves between-trials noise in sensory response to stimulus
su24  = unique(s(c==0.24));    %dir. for coh 24
su12  = unique(s(c==0.12));    %dir. for coh 24
su6   = unique(s(c==0.06));    %dir. for coh 24
Nsu24 = length(su24);
Nsu12 = length(su12);
Nsu6  = length(su6);
x1_24 = circ_vmrnd(su24(1)*ones(Ntrials,Ndots),fp(1)); %x for first direction
x1_12 = circ_vmrnd(su12(1)*ones(Ntrials,Ndots),fp(2));
x1_6  = circ_vmrnd(su6(1)*ones(Ntrials,Ndots),fp(3)); 

%calculate x for other directions s by rotating x so that new x densities peak at new s (speed up) 
%x new values are produced by rotating x values from s(1) by the angle between s(1) and s(new)
%calculate rotation from first direction for each direction
rotdeg24 = SLvectors2signedAngle(su24(1),su24,'radian');   %rotation angles that produce other s
rotdeg12 = SLvectors2signedAngle(su12(1),su12,'radian'); 
rotdeg6  = SLvectors2signedAngle(su6(1),su6,'radian');   
rotrad24 = SLde2r(rotdeg24,1);
rotrad12 = SLde2r(rotdeg12,1);
rotrad6  = SLde2r(rotdeg6,1);

%motor noise
motn = vmPdfs(0:1:359,0,km,'norm');

%get Pest each coherence
o24 = getPestEachCoh(Ntrials, Ndots, Ns, svec, su24, Nsu24 ,x1_24 ,rotrad24 ,0.24 ,fp(1),prior80,prior40,prior20,prior10,prand,motn);
o12 = getPestEachCoh(Ntrials, Ndots, Ns, svec, su12, Nsu12 ,x1_12 ,rotrad12 ,0.12 ,fp(2),prior80,prior40,prior20,prior10,prand,motn);
o6  = getPestEachCoh(Ntrials, Ndots, Ns, svec, su6 , Nsu6  ,x1_6  ,rotrad6  ,0.06 ,fp(3),prior80,prior40,prior20,prior10,prand,motn);

%Pestimates all conditions
%Pest is a Nestimates by Nconditions matrix
%[31*{24,80} 31*{24,40} 31*{24,20] 31*{24,10}  36*{12,80} .....]
%[Ndir * {coh * prior}      ......]
Pest = [o24.Pest o12.Pest o6.Pest];

%associated col conditions
sPest   = [o24.s_Pest   o12.s_Pest   o6.s_Pest  ];  %stim. directions
prPest  = [o24.pr_Pest  o12.pr_Pest  o6.pr_Pest ];  %prior 
cohPest = [o24.coh_Pest o12.coh_Pest o6.coh_Pest];  %coh
[mcond,posEx] = sortrows([prPest' cohPest' round(SLra2d(sPest))'],[1 2]);          %cond matrix

%-------------------------------
%Match data and model conditions
%-------------------------------
%get rid of non-existing conditions
[a,posNex] = setdiff(mcond,Cond,'rows');  %non existing conditions
modCond = mcond;
modCond(posNex,:) = [];                 
posEx(posNex)     = [];
Pest = Pest(:,posEx)';       %keep useful est densities

%sanity check
if size(modCond,1)~=size(Cond,1) || ~isempty(setdiff(modCond,Cond,'rows'))
    fprintf('Something"s wrong. The conditions associated with the model and the data do not match')
    keyboard
end

%get logl
for Ci = 1 : nCon

    Pest(Ci,Pest(Ci,:) <= 0) = 10^-320;         %get rid of <0 due to convolution
    Pest(Ci,:) = Pest(Ci,:)./sum(Pest(Ci,:),2); %proba
    logl(Ci) = sum(log(Pest(Ci,data(posC==Ci))));

end

loglAll = sum(logl);
fitP    = fp;


%get Probability of estimate for an individual coherence
function o = getPestEachCoh(Ntrials,Ndots,Ns,svec, su, Nsu,x1,rotrad,c,kappa,prior80,prior40,prior20,prior10,prand,motn)

%estimate space
predEst = 1 : 1 : 360;

%preallocate Probability of estimate given model
%each unique direction in coh 24%
x        = nan(Ntrials,Ndots,Nsu);    
coherent = nan(Ntrials,Ndots,Nsu); 
like     = nan(Ns,Ntrials);  
Pest_80 = nan(numel(predEst),Nsu); %estimate density
Pest_40 = nan(numel(predEst),Nsu);
Pest_20 = nan(numel(predEst),Nsu);
Pest_10 = nan(numel(predEst),Nsu);

%each unique stim dir this coh
parfor i = 1 : Nsu
    
    %x measurements   
    %to check that rotations are correct run this loop. ro should match rotdeg
    %repeat by changing the 3rd dim like:  x(:,i,2) to x(:,i,3).
    %x1 should match rotdeg(1)
    %x when i = 1 : Nsu is 2 should match rotdeg(2)
    %and so on ....
    %for j = 1 : size(x,2) 
    %   ro(:,j) = SLvectors2signedAngle(SLra2d(x(:,j),SLra2d(x1,'polar')
    %end
    x = slRotateMatrixValues(x1,rotrad(i),'dispOff');

    %coherent
    coherent = (rand(Ntrials,Ndots) < c);                       % which dots are moving coherently
    x = coherent .* x + (1-coherent) .* rand(Ntrials,Ndots)*2*pi;  % Ntrials by Ndots

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
    like_per_dot = (1 - c)/2/pi + c/2/pi/besseli(0,kappa,1) * exp(kappa * cos(bsxfun(@minus,x,svec))-kappa); %(1-c)/2pi + c*V(svec,x,kappa)
    like = prod(like_per_dot,3);

    %--------------
    %prior 80
    %--------------
    %0.01 s faster than bsxfun
    posterior_80 = like.*prior80(:,ones(1,Ntrials));
    
    %Mean
    xtemp_80 = sum(bsxfun(@times, cos(svec), posterior_80));
    ytemp_80 = sum(bsxfun(@times, sin(svec), posterior_80));
    postmeans_80 = atan2(ytemp_80,xtemp_80); % 1 by Ntrials)
    postmeans_80 = SLra2d(postmeans_80);

    %estimate distributions
    Pest_80(:,i) = hist(postmeans_80,predEst);
    Pest_80(:,i) = Pest_80(:,i)./sum(Pest_80(:,i));       %proba
    Pest_80(:,i) = (1 - prand)*Pest_80(:,i) + prand*(1/360);  %add random estimation
    Pest_80(:,i) = SLcircConv(Pest_80(:,i),motn)';           %add motor noise

    %---------
    %prior 40
    %---------
    %0.01 s faster than bsxfun
    posterior_40 = like.*prior40(:,ones(1,Ntrials));
    
    %Mean
    xtemp_40     = sum(bsxfun(@times, cos(svec), posterior_40));
    ytemp_40     = sum(bsxfun(@times, sin(svec), posterior_40));
    postmeans_40 = atan2(ytemp_40,xtemp_40); % 1 by Ntrials)
    postmeans_40 = SLra2d(postmeans_40);

    %estimate distributions
    Pest_40(:,i) = hist(postmeans_40,predEst);
    Pest_40(:,i) = Pest_40(:,i)./sum(Pest_40(:,i));       %proba
    Pest_40(:,i) = (1 - prand)*Pest_40(:,i) + prand*(1/360);  %add random estimation
    Pest_40(:,i) = SLcircConv(Pest_40(:,i),motn)';           %add motor noise

    %--------------
    %prior 20
    %--------------
    %0.01 s faster than bsxfun
    posterior_20 = like.*prior20(:,ones(1,Ntrials));    
    xtemp_20     = sum(bsxfun(@times, cos(svec), posterior_20));
    ytemp_20     = sum(bsxfun(@times, sin(svec), posterior_20));
    postmeans_20 = atan2(ytemp_20,xtemp_20);          %Mean
    postmeans_20 = SLra2d(postmeans_20);

    %estimate distributions
    Pest_20(:,i) = hist(postmeans_20,predEst);
    Pest_20(:,i) = Pest_20(:,i)./sum(Pest_20(:,i));       %proba
    Pest_20(:,i) = (1 - prand)*Pest_20(:,i) + prand*(1/360);  %add random estimation
    Pest_20(:,i) = SLcircConv(Pest_20(:,i),motn)';           %add motor noise


    %--------------
    %prior 10
    %--------------
    %0.01 s faster than bsxfun
    posterior_10 = like.*prior10(:,ones(1,Ntrials));
    
    %Mean
    xtemp_10 = sum(bsxfun(@times, cos(svec), posterior_10));
    ytemp_10 = sum(bsxfun(@times, sin(svec), posterior_10));
    postmeans_10 = atan2(ytemp_10,xtemp_10); % 1 by Ntrials)
    postmeans_10 = SLra2d(postmeans_10);

    %estimate distributions
    Pest_10(:,i) = hist(postmeans_10,predEst);
    Pest_10(:,i) = Pest_10(:,i)./sum(Pest_10(:,i));       %proba
    Pest_10(:,i) = (1 - prand)*Pest_10(:,i) + prand*(1/360);  %add random estimation
    Pest_10(:,i) = SLcircConv(Pest_10(:,i),motn)';           %add motor noise

end

%Pestimate
o.Pest     = [Pest_80 Pest_40 Pest_20 Pest_10];
o.s_Pest   = repmat(su',1,4);                 %each direction (col)
o.pr_Pest  = repmat([80 40 20 10],Nsu,1);     %each prior (col)
o.pr_Pest  = o.pr_Pest(:)';
o.coh_Pest = repmat(c,1,4*Nsu);               %each col coh