
%deriveLLHmargCoh.m
%
%     author: Steeve Laquitaine
%    purpose: derive LLH marginalized over coherences
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
%               o = SLderiveLLHmargCoh(180,1)
%


function output = SLderiveLLHmargCoh(Strue,nDots)

set(gcf,'color','w'); box off

%--------------------
%initialize variables
%--------------------
%an instance of dots directions in a 6% coh motion stimulus
%containing nDots (e.g.,100 dots). 6% coherent dots move in direction Strue
%(e.g., 225 deg) and 74% remaining dots in random directions.
px_SC = SLdeltaPdf(1:1:360,Strue,6,100);
x = randsample(1:1:360,nDots,true,px_SC);

%coherence
C = [0.06 ; 0.12 ; 0.24];   

%---------------------------------------------------------------
%scale to bypass numerical limits on prod of small probabilities
%---------------------------------------------------------------
%prevent prob numerical limit here due to prod small numbers!!!
pmax = 100;

%-----------------------
%Calculate Likelihood(S)
%-----------------------
for S = 1:360
    
    %P(C)
    nC = length(C);
    pC = repmat(1/nC,nC,1);
    
    %p(each x|S,C)
    nx = length(x);
    px_SC = nan(nx,nC);
    for i = 1 : nC
        
        %px_SC(:,i) = SLdeltaPdf(x,S,C(i)*pmax,pmax);
        
        vm = vmPdfs(1:1:360,S,10,'norm');
        px_SC(:,i) = (pmax - C(i))/360 + C(i)*pmax*vm(x);
        
    end
    
    %p(vector x|S,C)
    pVectx_SC = nan(nC,1);
    for i = 1 : nC
        pVectx_SC(i,:) = prod(px_SC(i,:),2);
    end
    
    %---------------------
    %marginalize coherence
    %---------------------
    %L(S) = p(vector x|S)...the likelihood of stim direction S.
    LS = sum( pVectx_SC .* pC, 1);
    
    fprintf(['(SLderiveLLH) LS = ',num2str(LS),'\n'])
    
    %plot
    hold all
    plot(S,LS,'k.')
    drawnow
end
xlabel('Stimulus direction space S (deg)')
ylabel('Likelihood(S)')
xlim([0 360])
text(180,log(LS)+0.1,'LS marginalized over the 3 coherences')

%outputs
output.dotsdirX = x;
output.coherences = C;
output.Strue = Strue;
output.nDots = nDots;

fprintf('(SLderiveLLHmargCoh)')