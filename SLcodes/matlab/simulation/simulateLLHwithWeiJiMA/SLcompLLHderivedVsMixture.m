

%SLcompLLHderivedVsMixture.m
%
%     author: Steeve Laquitaine
%    purpose: compare derived and mixture LLH for each coherence
%description:
%
%         inputs:
%               'Strue': True motion direction of stimulus
%               'nDots': number of dots in motion stimulus
%
%        outputs:
%            o.dotsdirx: simulated directions for each dot of the stimulus
%          o.coherences: 3 motion coherences which are marginalized over (6% 12% 24%)
%
%
%
%      usage:
%
%               o = SLcompLLHderivedVsMixture(75,1,[.8 .15 .05],'fixedEvidence',55)
%               o = SLcompLLHderivedVsMixture(75,1,[.33 .33 .33],'fixedEvidence',55)
%
%

function output = SLcompLLHderivedVsMixture(Strue,nDots,pC,varargin)

set(gcf,'color','w'); box off

%--------------------
%initialize variables
%--------------------
%if we sample direction evidence each time
if sum(strcmp(varargin,'sampledEvidence'))==1
    
    %an instance of dots directions in a 6% coh motion stimulus
    %containing nDots (e.g.,100 dots). 6% coherent dots move in direction Strue
    %(e.g., 225 deg) and 74% remaining dots in random directions.
    px_SC = SLdeltaPdf(1:1:360,Strue,6,100);
    x = randsample(1:1:360,nDots,true,px_SC);
end

%if we fix motion evidence
if sum(strcmp(varargin,'fixedEvidence'))==1
    pos = find(strcmp(varargin,'fixedEvidence')==1);
    x = varargin{pos+1};
end


%---------------------------------------------------------------
%scale to bypass numerical limits on prod of small probabilities
%---------------------------------------------------------------
%prevent prob numerical limit here due to prod small numbers!!!
pmax = 100;

%-----------------------
%Calculate Likelihood(S)
%-----------------------
Ci = [];
id = 0;
% rclr = [0 0 0; 0.5 0.5 0.5; 0.75 0.75 0.75];
clr = [1 0 0; 0.7 0.5 0.5; 0.5 0.5 0.5];

for C = 0.24; %[0.06 0.12 0.24]
    
    id = id + 1;
    
    for S = 1:360
        
        %p(each x|S,C)
        
        %px_SC(:,i) = SLdeltaPdf(x,S,C(i)*pmax,pmax);
        vm = vmPdfs(1:1:360,S,4,'norm');
        px_SC = (pmax - C)/360 + C*pmax*vm(x);
        
        %p(vector x|S,C)
        pVectx_SC = prod(px_SC,2);
        
        %L(S) = p(vector x|S,C)...the likelihood of stim direction S.
        for i = 1 : length(pC)
            LS(S) = pVectx_SC .* pC(i);
        end
        fprintf(['(SLderiveLLH) LS = ',num2str(LS),'\n'])
        
    end
    
    %plot
    subplot(1,2,1)
    hold all
    %1) LS = LS./sum(LS);
    plot(1:360,LS,'.','color',clr(id,:))
    xlabel('Stimulus direction space S (deg)')
    ylabel('Likelihood(S)')
    xlim([0 360])
    text(180,max(LS),[num2str(C),' % L(S) - coh'],'color',clr(id,:))
    title('Note that LLH do not sum to 1')
    
    %store
    Ci = [Ci;C];
    
    
    %combine with a prior for Bayesian inference
    %-------------------------------------------
    subplot(1,2,2)
    hold all
    prior = vmPdfs(1:1:360,225,0.7,'norm');
    post = (prior'.*LS)/(sum(prior'.*LS));
    plot(1:1:360,prior,'.b')
    plot(1:360,post,'.','color',clr(id,:))
    text(180,max(post)./2,[num2str(C),' % post - coh'],'color',clr(id,:))
    title('Posterior (sum to 1)')
    drawnow
    
end

legend

%outputs
output.dotsdirX = x;
output.coherences = Ci;
output.Strue = Strue;
output.nDots = nDots;



%% Why does mixture of von mises and uniform fail to create the two peaks
%with circular prior

%find the mixture LLH that matches best the derived LLH and test if it
%creates or not two peaks.
%My hypothesis is that the shape of the LLH are different.

%llh
kl=4;
kp=0.7;
wT=0.5;

llh = vmPdfs(1:1:360,x,kl,'norm');
uniform = (1/360)*ones(360,1);
llh = (1 - wT)*llh + wT*uniform;

subplot(1,2,1)
plot(1:1:360,llh,'color',[1 0.5 0.5])

%prior
prior = vmPdfs(1:1:360,225,kp,'norm');
subplot(1,2,2)
plot(1:1:360,prior,'--','color',[0.9 0.9 1])

%posterior
post = prior.*llh;
post = post/sum(post);
plot(1:1:360,post,'--','color',[1 0.5 0.5])


%%
figure;
hold all
plot(1:1:360,llh,'r')
plot(1:1:360,prior,'b')
plot(1:1:360,post,'k')



%%
ylim([0 0.06])
box off
hold all

%llh
llh = vmPdfs(dir,evidence(i),kl,'norm');
uniform = (1/360)*ones(numel(dir),1);
llh = (1 - weightTail)*llh  + weightTail*uniform;
area(llh,'facecolor',[.65 .65 .65],'linestyle','none')

%prior
prior = vmPdfs(dir,225,kp,'norm');
area(prior,'facecolor','b','linestyle','none')

%posterior
post = prior.*llh;
post = post/sum(post);
area(post,'facecolor','r','linestyle','none')

%std posterior
output.stdpost(i) = SLmakeStd(dir,post);

%legend
if i==1;
    legend('llh','prior','posterior')
    legend('boxoff')
end
ylim([0 max([llh(:);prior(:);post(:)])])
xlim([0 360])








