%author: steeve laquitaine
%date: 140424
%usage: NNwithXGain(10,16)
%for i=1:360;NNwithXGain(i); drawnow; end
%"motionDir" is the true motion direction. e.g., 10 degrees.
%"numNeurMTkj" is number of neurons in MT e.g., 16

%purpose: neural population simulation of different hypothetical effects of
%a prior on the population response.

%Description:
%Hardware
%- 16 MT neurons with equally spaced direction preferences
%- 45 degrees for a full width at half height of 106 degrees (Albright, 1984)
%For a Von Mises, it is a concentration parameter k of 2.34.
%- Encoding neurons are assumed statistically independent (Jazayeri, 2006)
%- Prior's std is k=2.75 (40 degrees std);
%- We assume same number of neurons and MT. Each neuron from MT
%projects to only one neuron in LIP.
%- decoding: 360 LIP neurons have 1 readout direction each.
%note: in Jazayeri, readout rule have cosine profiles..

%360 000 neurons in each voxel
%https://cfn.upenn.edu/aguirre/wiki/public:neurons_in_a_voxel

%Operations
%- The tuning curves fj in the sensory neurons (encoding) stay unchanged.
%- The logPrior adds to the PreSynaptic inputs; added to (Habenschuss,2013;
%Jazayeri,2006).
%- The weight of each encoder-decoder neuron is set according to
%Habenschuss et al., 2013 (consistent with Jazayeri et al., 2006)
%W_ed=logf_e(dir_d)+const
%- Winner-takes-all competition by divisive normalization(Carandini,Heeger)
%a common inhibitory signal with constant total firing rate in LIP
%(homeostasy).

%Questions
%- How are MT and LIP connected? Are there tuning curves in LIP?
%How does LIP represent posterior (decision variable), like MT?
%- Is prior gain represented in a probabilistic population of neurons with
%tuning curves (Ma et al., 2006)? See Simoncelli, 2009. No directly in the
%encoder-decoder synaptic weight after learning.


%To do
%simulate MT and LIP voxels responses (10-20 voxels);


function [MTBOLDmnByPref,LIPBOLDmnByPref,k,motionspace,VoxPrefDirMT,VoxPrefDirLIP]=NNwithXGain(motionDir,...
    priorstd,...
    priormean,...
    numVoxelsMT,...
    numVoxelsLIP,...
    NeurInaVoxelMT,...
    NeurInaVoxelLIP,...
    numNeurMTkj,...
    OptiondisplayNetwork,...
    OptiondisplayBOLD)

%parameters
MTBOLDmn=[];
LIPBOLDmn=[];
priorstd=str2double(priorstd(find(priorstd=='=')+1:end));
numVoxelsMT=str2double(numVoxelsMT(find(numVoxelsMT=='=')+1:end));
numVoxelsLIP=str2double(numVoxelsLIP(find(numVoxelsLIP=='=')+1:end));
NeurInaVoxelMT=str2double(NeurInaVoxelMT(find(NeurInaVoxelMT=='=')+1:end));
NeurInaVoxelLIP=str2double(NeurInaVoxelLIP(find(NeurInaVoxelLIP=='=')+1:end));
motionspace=1:1:360;

%LIP decoder neurons
k=1:1:360;
numNeurLIPk=numel(k);

%MT's encoder (preferred directions must be evenly spaced !!)
DirPref=1:360/numNeurMTkj:360;
if fix(360/numNeurMTkj)-360/numNeurMTkj~=0
    fprintf('%30s \n','Number of MT neurons must be multiple of 360')
    keyboard
end
DirK=2.34;
Fkj=vmPdfs(motionspace,DirPref',DirK(ones(numel(DirPref),1),1),'norm');

%MT's output firing responses 
rMTj=nan(numel(motionDir),numNeurMTkj);

%MT's neurons in each voxel whichNeurInaVoxel (neurons,voxels)
whichNeurInaVoxel=nan(numel(NeurInaVoxelMT),numVoxelsMT);
for j=1:numVoxelsMT
    whichNeurInaVoxel(:,j)=randsample(1:1:numNeurMTkj,numel(NeurInaVoxelMT),true);
end
MTBOLD=nan(numVoxelsMT,1);

%LIP total output firing response (homeostasy)
rtotal=10000;

%LIP's neurons in each voxel whichNeurInaVoxelLIP (neurons,voxels)
whichNeurInaVoxelLIP=nan(numel(NeurInaVoxelLIP),numVoxelsLIP);
for j=1:numVoxelsLIP
    whichNeurInaVoxelLIP(:,j)=randsample(1:1:numNeurLIPk,numel(NeurInaVoxelLIP),true);
end
LIPBOLD=nan(numVoxelsLIP,1);

%Prior
%PriorMean=225;
PriorMean=str2double(priormean(find(priormean=='=')+1:end));
PriorGain=vmPdfs(1:1:360,PriorMean,priorstd,'norm');


%Simulate firing and BOLD responses for each motion direction
fprintf('%50s \n','Now running the network...if displayBOLD is on wait for upcoming plot')
for i=1:numel(motionDir);
    
    %print trial
    fprintf('%12s \n',num2str(i))
    
    %MT's firing responses
    dirThisTrial=k==motionDir(i);
    rMTj(i,:)=Fkj(dirThisTrial,:);
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        figure(1)
        screen=get(0,'ScreenSize');
        set(gcf,'position', [0*screen([3 4]) 0.4*screen(3) screen(4)], 'color','w')
        clf
        subplot(531);
        C=linspecer(numNeurMTkj);
        set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
        hold all
        plot(Fkj,'linesmoothing','on','linewidth',3)
        pos=motionspace(motionspace==motionDir(i));
        plot(pos(ones(1,numNeurMTkj)),rMTj(i,:),'.','markersize',20)
        xlim([0 360])
        xlabel('Motion direction (degrees)','fontsize',14)
        ylabel('Firing rate (sp/sec)','fontsize',14)
        title('MT direction tuned neurons','fontsize',14,'fontweight','bold')
        box off
    end
    
    
    %Prior with multiplicative gain
    %Need to think about how to relate Prior's std to the multiplicative gain
    %but for now multiplicative gain's std is prior's std.
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        sp(1)=subplot(539);
        plot(log(PriorGain),'color',[1 0 0],'linesmoothing','on','linewidth',3)
        xlim([0 360])
        xlabel('Motion direction (degrees)','fontsize',14)
        ylabel('Weight','fontsize',14)
        title({'Prior weight (hypothetical)','logprior(\theta_k)'},...
            'fontsize',14,'fontweight','bold','color','r')
        box off
    end
        
    %Synaptic weights of each MT neuron - LIP neurons. It is the
    %neurons logLlh weight (Jazayeri et al., 2006). Some terms isolated
    %during posterior derivation are assumed constant over motion directions
    %and have been dropped(see Jazayeri et al., 2006).
    Wkj=log(Fkj);
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        subplot(534)
        set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
        plot(k,Wkj,'linesmoothing','on','linewidth',3)
        xlim([0 360])
        title({'Connections weights','(logf_j)'},...
            'fontsize',14,'fontweight','bold')
        xlabel('Motion direction (degrees)','fontsize',14)
        ylabel('Weight','fontsize',14)
        box off
    end
   
    
    %Pooled synaptic weights (FR, Wkj, prior)
    %Synaptic weights are pooled for each readout neuron.
    %example
    %rMTj(i,:)=Fkj(k==motionDir(i),:);
    rMTjThisTrial=rMTj(i,:);
    rMTkj=rMTjThisTrial(ones(numel(k),1),:);
    pooledWj=sum(rMTkj.*Wkj,2);
    %pooledWj=sum(Wkj,2);
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        sp(2)=subplot(537);
        plot(k,pooledWj,'k','linesmoothing','on','linewidth',3)
        xlim([0 360])
        ylim([1.02*min(pooledWj) 0.98*max(pooledWj)])
        title({'Pooled presynaptic inputs','(logllh(\theta_k)=\Sigma_j r_j*logf_j(\theta_k))'},...
            'fontsize',14,'fontweight','bold')
        xlabel({'Readout neurons','by readout direction (degrees)'},'fontsize',14)
        ylabel('Weight','fontsize',14)
        box off
    end
    
    %LIP's decoder neurons output firing rate(max pooling)
    %Based on Habenschuss' derivation assuming divisive normalization
    %(inhibition). It is a softmax or attractor network. It makes sure the total
    %firing rate in LIP remains constant.
    %- if a readout neuron firing rate is r_k=exp(u_k)
    %- if r_k depends on u_k=sum_j(w_k.*r_j + logprior - I)
    %- if r_total=sum_k(r_k) is maintained
    %- the ideal divisive inhibition is I=log(Sum(exp(w_k.*r_j + logprior))) -
    %logr_total
    %note: pooledWj is w_k.*r_j.
    %r_j is neuron j's firing
    %and
    %in Habenschuss beta=1;
    beta=1;
    rLIPk=rtotal.*exp(beta*(pooledWj+log(PriorGain)))./sum(exp(beta*(pooledWj+log(PriorGain))));
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        sp(3)=subplot(5,3,11);
        plot(k,rLIPk,'k','linesmoothing','on','linewidth',3)
        xlim([0 360])
        ylim([0.75*min(rLIPk) 1.5*max(rLIPk)])
        title({'LIP neurons responses','logllh(\theta_k)+logprior'},...
            'fontsize',14,'fontweight','bold')
        xlabel({'Readout neurons','by readout direction (degrees)'},'fontsize',14)
        ylabel('Firing rate (sp/sec)','fontsize',14)
        box off
    end
    
    
    %look at prior and llh weights at the same scale
    l=log(PriorGain); 
    minY=min([l(:);pooledWj(:)]);
    maxY=max([l(:);pooledWj(:)]);
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        sp(4)=subplot(538);
        hold all
        plot(k,log(PriorGain),'color',[1 0 0],'linesmoothing','on','linewidth',3)
        plot(k,pooledWj,'color','k','linesmoothing','on','linewidth',3)
        xlim([0 360])
        xlabel('Motion direction (degrees)','fontsize',14)
        ylabel('Weight','fontsize',14)
        title({'logPrior and logllh at same scale'},...
            'fontsize',14,'fontweight','bold','color','k')
        box off
        set(sp(4),'ylim',[minY-0.1*(maxY-minY) maxY+0.1*(maxY-minY)])
    end

    %read-out direction
    [maxrLIPk,pos]=max(rLIPk);
    estimate=k(pos);
    if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
        subplot(5,3,14)
        drawVectors(motionDir(i),1)
        hold on
        drawVectors(estimate,1)
        hold on
        drawVectors(PriorMean,1)
        coor=polar2cartesian(estimate,2.1);
        text(coor(1),coor(2),num2str(estimate),'color','r','fontsize',20,...
            'fontweight','bold');
    end   
    
    % fMRI BOLD prediction for MT
    %MT only produces sensory responses r_j.
    %One MT voxel is made of 100 direction-tuned neurons randomly sampled from
    %the initial pool of 16 MT neurons.
    for j=1:numVoxelsMT
        MTresponsesInaVoxel=rMTjThisTrial(whichNeurInaVoxel(:,j));
        MTBOLD(j)=sum(MTresponsesInaVoxel);
    end
    
    % fMRI BOLD prediction for LIP
    %LIP produces percept responses rLIPk.
    %One LIP voxel is made of 1000 readout neurons randomly sampled from
    %the initial pool of 360 LIP neurons.
    for j=1:numVoxelsLIP
        LIPresponsesInaVoxel=rLIPk(whichNeurInaVoxelLIP(:,j));
        LIPBOLD(j)=sum(LIPresponsesInaVoxel);
    end
    MTBOLDmn=[MTBOLDmn MTBOLD];
    LIPBOLDmn=[LIPBOLDmn LIPBOLD];
    drawnow
end
    
%print message
fprintf('%50s \n','Wait...now sorting voxel by direction selectivity and preparing plot...')

%calculate voxels preferred direction
VoxPrefDirMT=nan(numVoxelsMT,1);

%progress indicator
h=waitbar(0,'Please wait...');
for vox=1:numVoxelsMT;
    [m,s]=makeStat(MTBOLDmn(vox,:)',motionDir);
    [~,PrefDirMT]=max(m);
    VoxPrefDirMT(vox)=PrefDirMT;
    %progress indicator
    perc=fix((vox/numVoxelsMT)*100)/100;
    waitbar(perc,h,sprintf('%d%% sorting MT voxels based on direction selectivity',perc*100))
end
close(h)

%progress indicator
h=waitbar(0,'Please wait...');
VoxPrefDirLIP=nan(numVoxelsLIP,1);
for vox=1:numVoxelsLIP;
    [m,s]=makeStat(LIPBOLDmn(vox,:)',motionDir);
    [~,PrefDirLIP]=max(m);
    VoxPrefDirLIP(vox)=PrefDirLIP;
    %progress indicator
    perc=fix((vox/numVoxelsLIP*100))/100;
    waitbar(perc,h,sprintf('%d%% sorting LIP voxels based on direction selectivity',perc*100))
end
close(h)


%draw MT and LIP BOLD responses sorted by voxels preferred directions for
%each motion directions displayed
%sort responses by voxels preferred directions
[VoxPrefDirMT,IMT]=sort(VoxPrefDirMT);
MTBOLDmnByPref=MTBOLDmn(IMT,:);

[VoxPrefDirLIP,ILIP]=sort(VoxPrefDirLIP);
LIPBOLDmnByPref=LIPBOLDmn(ILIP,:);

if strcmp(OptiondisplayBOLD,'displayBOLD=on')==1
    figure;
    screen=get(0,'ScreenSize');
    set(gcf,'position', [0.4*screen(3) 0.9*screen(4) 0.4*[screen(3) screen(4)]], 'color','w')
    unqDir=unique(motionDir);
    for i=1:numel(unqDir);
        
        %MT Bold
        subplot(211)
        meanMTBOLDoverRep=mean(MTBOLDmnByPref(:,motionDir==unqDir(i)),2);
        plot(1:1:numVoxelsMT,meanMTBOLDoverRep,'color',[.5 .5 .5],'linesmooth','on','linewidth',2)
        title({'MT','BOLD=\Sigma responses of sampled 100000 neurons'},'fontsize',14)
        xlabel('Voxels by preferred directions (degrees)','fontsize',14)
        ylabel('Simulated Bold response','fontsize',14)
        xlim([0 numVoxelsMT])
        ylim([min(MTBOLDmn(:)) max(MTBOLDmn(:))])
        set(gca,'fontsize',14,'xtick',1:10:numVoxelsMT,'xticklabel',VoxPrefDirMT(1:10:end))
        box off
        
        %LIP BOLD
        subplot(212)
        meanLIPBOLDoverRep=mean(LIPBOLDmnByPref(:,motionDir==unqDir(i)),2);
        plot(1:1:numVoxelsLIP,meanLIPBOLDoverRep,'color',[.5 .5 .5],'linesmooth','on','linewidth',2)
        title({'LIP','BOLD=\Sigma responses of sampled 100000 neurons'},'fontsize',14)
        xlabel('Voxels by preferred directions (degrees)','fontsize',14)
        ylabel('Simulated Bold response','fontsize',14)
        xlim([0 numVoxelsLIP])
        ylim([min(LIPBOLDmn(:)) max(LIPBOLDmn(:))])
        set(gca,'fontsize',14,'xtick',1:10:numVoxelsLIP,'xticklabel',VoxPrefDirLIP(1:10:end))
        box off
        drawnow
    end
end

%density of voxels preferred directions
if strcmp(OptiondisplayBOLD,'displayBOLD=on')==1
    figure('color','w')
    set(gcf,'position', [0.4*screen(3) 0.05*screen(4) 0.4*[screen(3) screen(4)]], 'color','w')
    subplot(121)
    title('MT voxels')
    hist(VoxPrefDirMT)
    set(get(gca,'child'),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
    box off
    xlabel('Voxels preferred directions (degrees)')
    ylabel('Number of voxels')

    subplot(122)
    title('LIP voxels')
    hist(VoxPrefDirLIP)
    set(get(gca,'child'),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
    xlabel('Voxels preferred directions (degrees)')
    ylabel('Number of voxels')
    box off
end


