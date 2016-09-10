%NNwithXGain2.m

% author: steeve laquitaine
%   date: 140424 last modification 140602
%  
%  usage: 
%       dirTrained=repmat(1:1:360,1,1);
%       whichNeurInaVoxMT=[];
%       neuronIDinVoxPOST=[];
%       neuronIDinVoxPrior=[];
%       [llhBOLDmn,postBOLDmn,priorBold]=NNwithXGain2(dirTrained,...
%           'priorstd=0',...
%           'priormean=225',...
%           'numVoxsMT=1',...
%           'numVoxInPosteriorArea=100',...
%           'NeurInaVoxMT=1',...
%           'numNeuroInaVoxPOST=1000',...
%           whichNeurInaVoxMT,...
%           neuronIDinVoxPOST,...
%           neuronIDinVoxPrior,...
%           15,...
%           'sortVoxels=on',...
%           'displayNetwork=on',...
%           'displayBOLD=off');
% 
% purpose: Simulation of neural and BOLD responses of a Bayesian network of
%          llh, prior and posterior areas
% 
% Description:
%         Hardware
%         - 16 llh neurons with equally spaced direction preferences
%         - 45 degrees for a full width at half height of 106 degrees (Albright, 1984)
%         For a Von Mises, it is a concentration parameter k of 2.34.
%         - Encoding neurons are assumed statistically independent (Jazayeri, 2006)
%         - Prior's concentration is k=0.1;
%         - We assume same number of neurons and llh. Each neuron from llh
%         projects to only one neuron in Post.
%         - decoding: 360 Post neurons have 1 readout direction each.
%         note: in Jazayeri, readout rule have cosine profiles..
% 
%         360 000 neurons in each voxel
%         https://cfn.upenn.edu/aguirre/wiki/public:neurons_in_a_voxel
% 
% Operations
%         - The tuning curves fj in the sensory neurons (encoding) stay unchanged.
%         - The logPrior adds to the PreSynaptic inputs; added to (Habenschuss,2013;
%         Jazayeri,2006).
%         - The weight of each encoder-decoder neuron is set according to
%         Habenschuss et al., 2013 (consistent with Jazayeri et al., 2006)
%         W_ed=logf_e(dir_d)+const
%         - Winner-takes-all competition by divisive normalization(Carandini,Heeger)
%         a common inhibitory signal with constant total firing rate in Post
%         (homeostasy).
% 
% Questions
%         - How are llh and Post connected? Are there tuning curves in Post?
%         How does Post represent posterior (decision variable), like llh?
%         - Is prior gain represented in a probabilistic population of neurons with
%         tuning curves (Ma et al., 2006)? See Simoncelli, 2009. No directly in the
%         encoder-decoder synaptic weight after learning.
%
% To do:
% organize the code with first, neural network responses and second, BOLD responses

function [llhBOLDmn,PostBOLDmn,PriorBOLDmn,readoutTunings,...
    VoxRespAmpPrior,VoxPrefDirllh,VoxPrefDirPost,neuronIDinVoxLLH,...
    neuronIDinVoxPOST,neuronIDinVoxPrior,output] = NNwithXGain2(motionDir,...
    priorstd,...
    priormean,...
    numVoxInLLHarea,...
    numVoxInPosteriorArea,...
    numNeuroInaVoxLLH,...
    numNeuroInaVoxPOST,...
    neuronIDinVoxLLH,...
    neuronIDinVoxPOST,...
    neuronIDinVoxPrior,...
    numTuningInLLH,...
    sortVoxels,...
    OptiondisplayNetwork,...
    OptiondisplayBOLD)


%-------------------------------------------
%Network architecture and neurons' properties
%-------------------------------------------

%Readout neurons (posterior) 
%---------------------------
%readout neurons
output.readoutTunings = 1:1:360;
numReadoutNeur = numel(output.readoutTunings);


%sensory neurons (evidence)
%--------------------------
%neurons' tuning in sensory evidence-area 
%note: preferred directions must be evenly spaced
%direction selectivity
output.sensTunings = 1:360/numTuningInLLH:360;

%warning
if fix(360/numTuningInLLH) - 360/numTuningInLLH~=0
    fprintf('%s \n','(NNwithXGain2) Number of llh neurons must be multiple of 360')
    keyboard
end

%tuning width and shape
DirK = 2.34;
output.sensNeurTuning = vmPdfs(output.readoutTunings,output.sensTunings',DirK(ones(numel(output.sensTunings),1),1),'norm');

%spike counts scaled to 200 spikes max
%working in spike counts (integer is important for following operations)
%add lapse rate of 1 spike because we use logarithm later
output.sensNeurTuning = 16*output.sensNeurTuning./max(output.sensNeurTuning(:));
output.sensNeurTuning = round(output.sensNeurTuning*10) + 1;
output.sensNeurTuning = output.sensNeurTuning';

%sensory presynaptic inputs 
%basal sensory presynaptic input
%(Jazayeri et al., 2006)
fprintf('%s \n','Sensory neurons send basal presyaptic inputs to readout neurons')
output.weightSensPreSynTuning = log(output.sensNeurTuning);
display5(OptiondisplayNetwork,output)


%prior neurons (prior)
%---------------------
%neurons' basic properties
%This prior representation model is same as sensory evidence but with constant input
%motion direction at prior mode (not very intuitive though).
%A more intuitive scenario is one where the prior is implemented 
%in readout neurons' synapses (e.g., in LIP).
output.PriorMean = str2double(priormean(find(priormean=='=')+1:end));

%tuning in spike counts 
%scaled to 160 spikes max because 170 is the nax factorial function can
%handle (we use it later for normalization)
output.priorNeurTuning = vmPdfs(output.readoutTunings,output.sensTunings',DirK(ones(numel(output.sensTunings),1),1),'norm');
output.priorNeurTuning = 16*output.priorNeurTuning./max(output.priorNeurTuning(:));
output.priorNeurTuning = round(output.priorNeurTuning*10);

%strength
%strength is a multiplicative gain. It does not change neurons' tuning widths
%When prior is most certain (std = 0), tuning curves amplitude is twice
%higher. On the opposite weakest prior (strength = 0) have flat tuning curves.
output.PriorStd = str2double(priorstd(find(priorstd=='=')+1:end));
strength = 2./(output.PriorStd.^2+1);
output.priorNeurTuning = strength.*output.priorNeurTuning';

%scale tuning to bypass computational limits
%trick: adding a lapse firing rate of 1 spike permits to set prior weight at 0
%when prior is flat and neurons responses at 
%trick: We normalize max tuning responses to 171 spikes. 320 spikes is the max
%spikes we get when prior is strongest (std=0) and tuning curves are
%narrowest. We normalize tuning response such that tuning curves are twice 
%higher for this condition than in basal condition (strength=1);
output.priorNeurTuning = 170*output.priorNeurTuning/320;
output.priorNeurTuning = round(output.priorNeurTuning) + 1;

%responses (spikes)
%(firing "as if" motion direction was always the prior mean)
%note: neurons selective for prior more respond more
output.priorNeurFIRING = output.priorNeurTuning(:,output.readoutTunings==output.PriorMean);

%Prior's presynaptic weights (built-in)
output.weightPRIORPreSynTuning = log(output.priorNeurTuning);


%prior synaptic weights at readout neurons
%-----------------------------------------
%They come from pooling together:
%   - sensory neurons' current responses to the motion
%   - sensory neurons' built-in tuning weights (output.weightSensPreSynTuning)
%   - prior weights
%
%Synaptic weights are pooled for each readout neuron.
%note: The max weight should always peak for readout neurons selective to prior mode 
%
%The additional factors are normalization factors to make sure weights are
%positive and in a correct range. Those are very important if we consider modulating 
%prior strength!
%
%note max variable for factorial is 170
%This is Jazayeri et al. 2006's exact equation
%Warning;
if max(output.priorNeurFIRING)>170
    fprintf('%s \n',['Prior neurons responses are too high (>170)',...
        ' factorial can not keep up with it. Please SCALE DOWN !!'])
    keyboard
end

fprintf('%s \n','Calculating Prior presynaptic inputs....')
for readoutn = 1 : numReadoutNeur
    output.weightPRIOR(readoutn) = sum(output.priorNeurFIRING.*output.weightPRIORPreSynTuning(:,readoutn),1) ...
        - sum(output.weightPRIORPreSynTuning(:,readoutn)) ...
        - sum(log(factorial(output.priorNeurFIRING)));
end
output.weightPRIOR = output.weightPRIOR';



%------------------------
%ASSIGN NEURONS TO VOXELS
%------------------------
%sensory evidence
%----------------
%llh's neurons in each voxel neuronIDinVoxLLH (neurons,voxels). You can
%input the neurons yourself e.g., when you want to run different 
%analyses on the same voxels' set. If you sample from a pool of 15 direction-selective 
%neurons, neuronIDinVoxLLH values range between 1:1:15.
%If you dont' input any value the code automatically simulate a set of
%voxels
numNeuroInaVoxLLH = str2double(numNeuroInaVoxLLH(find(numNeuroInaVoxLLH=='=')+1:end));
numVoxInLLHarea = str2double(numVoxInLLHarea(find(numVoxInLLHarea=='=')+1:end));
if isempty(neuronIDinVoxLLH)
    neuronIDinVoxLLH = nan(numNeuroInaVoxLLH,numVoxInLLHarea);
    for j = 1 : numVoxInLLHarea
        neuronIDinVoxLLH(:,j) = randsample(1:1:numTuningInLLH,numNeuroInaVoxLLH,true);
    end
else
    fprintf('%s \n','(NNwithXGain2) Neurons have been assigned to their llh voxels..done')
end
output.motionDir = motionDir;
llhBOLDmn = nan(numVoxInLLHarea,numel(output.motionDir));

%Posterior
%---------
%Same idea as llh with 360 neurons.
numNeuroInaVoxPOST = str2double(numNeuroInaVoxPOST(find(numNeuroInaVoxPOST=='=')+1:end));
numVoxInPosteriorArea = str2double(numVoxInPosteriorArea(find(numVoxInPosteriorArea=='=')+1:end));
if isempty(neuronIDinVoxPOST)
    neuronIDinVoxPOST = nan(numNeuroInaVoxPOST,numVoxInPosteriorArea);
    for j=1:numVoxInPosteriorArea
        neuronIDinVoxPOST(:,j) = randsample(1:1:numReadoutNeur,numNeuroInaVoxPOST,true);
    end
else
    fprintf('%s \n','(NNwithXGain2) Neurons have been assigned to their posterior voxels...done')
end
PostBOLDmn = nan(numVoxInPosteriorArea,numel(output.motionDir));

%prior
%-----
numVoxsPrior = 1;
numNeuroInVoxPRIOR = numel(neuronIDinVoxPrior);
if isempty(neuronIDinVoxPrior)
    neuronIDinVoxPrior=nan(numNeuroInVoxPRIOR,numVoxsPrior);
    for j=1:numVoxsPrior
        neuronIDinVoxPrior(:,j)=randsample(1:1:numReadoutNeur,numNeuroInVoxPRIOR,true);
    end
else
    fprintf('%s \n','(NNwithXGain2) Neurons have been assigned to their prior voxels...done')
end
PriorBOLDmn = nan(numVoxsPrior,numel(output.motionDir));


%---------------
%Run the network
%---------------
%Initialize
fprintf('%s \n','(NNwithXGain2) Now running the Bayesian network...if "displayBOLD==on" wait for plot')
fprintf('%s \n','Likelihood inputs and readouts are updated at the display of each motion direction....')
screen = get(0,'ScreenSize');

%Graphics
C = linspecer(numTuningInLLH);

%neurons' responses
output.sensNeurFIRING = nan(numTuningInLLH,numel(output.motionDir));

%motion directions
for i = 1 : numel(output.motionDir);
    
    %progress
    SLprintProgress(i,numel(output.motionDir))
    
    %motion direction
    fprintf('%s \n',['Display direction ',num2str(output.motionDir(i)),' deg'])
    
    %sensory neurons response
    %------------------------
    fprintf('%s \n','Sensory neurons respond to motion (spikes)')
    dirThisTrial = output.readoutTunings==output.motionDir(i);
    output.sensNeurFIRING(:,i) = output.sensNeurTuning(:,dirThisTrial); 
    display3(OptiondisplayNetwork,output,numTuningInLLH,i)
    
    %likelihood presynaptic weights 
    %sensory presynaptic input after sensory neurons' response to stimulus
    %They come from pooling together:
    %   - sensory neurons' current responses to the motion with
    %   - sensory neurons' basal sensory presynaptic input and
    %   - prior presynaptic inputs
    %Synaptic inputs are pooled at each readout neuron.
    %This is Jazayeri et al. 2006's exact equation
    for readoutn = 1 : numReadoutNeur
        output.weightLLH(readoutn,i) = sum(output.sensNeurFIRING(:,i).*output.weightSensPreSynTuning(:,readoutn),1) ...
           - sum(output.weightSensPreSynTuning(:,readoutn)) ...
           - sum(log(factorial(output.sensNeurFIRING(:,i))));
    end    
    display6(OptiondisplayNetwork,output,i)

    %Prior 
    %-----
    display4(OptiondisplayNetwork,output)
    
    %prior and llh inputs on same scale
    display8(OptiondisplayNetwork,output,i)
    
    %Posterior
    %---------
    %1 - integrate prior and llh weights
    %2 - divisive competition between readout neurons
    %
    %Post's decoder neurons output firing rate(max pooling)
    %Based on Habenschuss' derivation assuming divisive normalization
    %(inhibition). It is a softmax or attractor network. It makes sure the total
    %firing rate in Post remains constant.
    %- if a readout neuron firing rate is r_k=exp(u_k)
    %- if r_k depends on u_k=sum_j(w_k.*r_j + logprior - I)
    %- if r_total=sum_k(r_k) is maintained
    %- the ideal divisive inhibition is I=log(Sum(exp(w_k.*r_j + logprior))) -
    %logr_total
    %note: weightLLH(:,i) is w_k.*r_j.
    %r_j is neuron j's firing
    %and
    %in Habenschuss beta=1;
    
    %readout neurons' responses in spikes
    %exponential is good because response cannot be negative
    %we scale with absolute value becuase scaling shouldn't change the
    %sign. "beta" must be relatively low otherwise we reach computational
    %limits with exponential. "rtotal" is large to have visible enough response
    %amplitude.
    beta = 0.01;
    rtotal = 10000;
    output.rPostk(:,i) = rtotal.*exp(beta*(output.weightLLH(:,i) + output.weightPRIOR))./...
        abs(sum(exp(beta*(output.weightLLH(:,i) + output.weightPRIOR))));
    output.rPostk(:,i) = round(output.rPostk(:,i));
    display7(OptiondisplayNetwork,output,i)
 
    %read-out the direction
    [~,pos] = max(output.rPostk(:,i));
    output.estimate(i) = output.readoutTunings(pos);
    display9(OptiondisplayNetwork,output,i)
   
    %--------------------
    %make Voxels and BOLD
    %--------------------
    %llh
    %llh only produces sensory responses r_j. A voxel is made out of e.g., 1000 
    %neurons randomly sampled from the initial pool of 15 llh direction-
    %selective neurons.
    llhBOLD = makeVoxels(neuronIDinVoxLLH,output.sensNeurFIRING(:,i),'verbose=off');    
    llhBOLDmn(:,i) = llhBOLD;
    
    %Posterior
    %Post produces responses output.rPostk. An Post voxel is made out of 1000 readout 
    %neurons randomly sampled from the initial pool of 360 Post neurons
    %selective for one direction only.
    PostBOLD = makeVoxels(neuronIDinVoxPOST,output.rPostk(:,i),'verbose=off');
    PostBOLDmn(:,i) = PostBOLD;
    
    %Prior
    %We expect BOLD to be same across displayed directions because the prior gain
    %does not depend on displayed motion directions. The firing responses
    %of the sampled neurons in the area is same across presentation of
    %motion direction and only reflect prior gain.
    PriorBOLD = makeVoxels(neuronIDinVoxPrior,output.priorNeurFIRING,'verbose=off');
    PriorBOLDmn(:,i) = PriorBOLD;
    drawnow
end


%end progress
fprintf('\n');


%VOXELS
VoxPrefDirllh=[];
VoxPrefDirPost=[];
VoxRespAmpPrior=[];

%case sort voxel (TO BE DEBUGGED)
if strcmp(sortVoxels,'sortVoxels=on')==1

    %print
    fprintf('%s',['(NNwithXGain2) Getting voxels direction',...
        ' selectivity and sorting voxels by direction selectivity...'])

    %Calculate average BOLD for each direction and voxels preferred
    %direction for llh and Post. BOLD is averaged over repetition of directions.
    %Preferred direction is the direction that elicits maximum average BOLD.
    %h=waitbar(0,'Please wait...'); Prior voxels have no direction
    %selectivity and fire the same way whatever displayed direction.
    VoxPrefDirllh = nan(numVoxInLLHarea,1);
    %     parfor vox=1:numVoxInLLHarea;
    for vox=1:numVoxInLLHarea;
        
        [m,s]=makeStat(llhBOLDmn(vox,:)',output.motionDir);
        [~,PrefDirllh]=max(m);
        VoxPrefDirllh(vox)=PrefDirllh;
        
        %progress
        %perc=fix((vox/numVoxInLLHarea)*100)/100;
        %waitbar(perc,h,sprintf('%d%% sorting llh voxels based on direction selectivity',perc*100))
    end
    %close(h)
    
    %progress
    %h=waitbar(0,'Please wait...');
    VoxPrefDirPost=nan(numVoxInPosteriorArea,1);
    %     parfor vox=1:numVoxInPosteriorArea;
    for vox=1:numVoxInPosteriorArea;
        
        [m,s]=makeStat(PostBOLDmn(vox,:)',output.motionDir);
        [~,PrefDirPost]=max(m);
        VoxPrefDirPost(vox)=PrefDirPost;
        
        %progress
        %perc=fix((vox/numVoxInPosteriorArea*100))/100;
        %waitbar(perc,h,sprintf('%d%% sorting Post voxels based on direction selectivity',perc*100))
    end
    %close(h)
    
    %sort matrix of BOLD by voxels preferred directions
    [VoxPrefDirllh,Illh] = sort(VoxPrefDirllh);
    llhBOLDmn = llhBOLDmn(Illh,:);
    
    [VoxPrefDirPost,IPost] = sort(VoxPrefDirPost);
    PostBOLDmn = PostBOLDmn(IPost,:);
    
    [VoxRespAmpPrior,IPrior] = sort(PriorBOLDmn(:,1));
    PriorBOLDmn = PriorBOLDmn(IPrior,:);
   
    %stick neurons to their voxels
    neuronIDinVoxLLH = neuronIDinVoxLLH(:,Illh);
    neuronIDinVoxPOST = neuronIDinVoxPOST(:,IPost);
    neuronIDinVoxPrior = neuronIDinVoxPrior(:,IPrior);
    fprintf('%s \n','done')
end

%draw llh and Post BOLD.
if strcmp(OptiondisplayBOLD,'displayBOLD=on')==1
    
    %calculate average BOLD response for unique direction
    unqDir=unique(output.motionDir);
    for i=1:numel(unqDir);
        %llh, posterior and Prior
        meanllhBOLDoverRep(:,i)=mean(llhBOLDmn(:,output.motionDir==unqDir(i)),2);
        meanPostBOLDoverRep(:,i)=mean(PostBOLDmn(:,output.motionDir==unqDir(i)),2);
        meanPriorBOLDoverRep(:,i)=mean(PriorBOLDmn(:,output.motionDir==unqDir(i)),2);
    end
    
    %draw BOLD responses for each voxels for unique directions (each screen
    %refresh)
    figure;
    set(gcf,'position',...
        [0.4*screen(3) 0.9*screen(4) 0.4*[screen(3) screen(4)]],...
        'color','w')
    for i=1:numel(unqDir);
        
        %llh
        subplot(311)
        plot(1:1:numVoxInLLHarea,meanllhBOLDoverRep(:,i),'color',[.5 .5 .5],...
            'linesmooth','on','linewidth',2)
        title({'llh',['BOLD=\Sigma FR of',num2str(numNeuroInaVoxLLH),...
            ' sampled neurons for direction ',num2str(unqDir(i)),' degrees']},...
            'fontsize',14)
        xlabel('Voxels by preferred directions (degrees)','fontsize',14)
        ylabel('Simulated Bold response','fontsize',14)
        xlim([0 numVoxInLLHarea])
        ylim([0.99*min(meanllhBOLDoverRep(:)) 1.01*max(meanllhBOLDoverRep(:))])
        
        if strcmp(sortVoxels,'sortVoxels=on')==1 
            set(gca,'fontsize',14,'xtick',1:10:numVoxInLLHarea,'xticklabel',...
                VoxPrefDirllh(1:10:end))
        end
        box off
        
        %Post
        subplot(312)
        plot(1:1:numVoxInPosteriorArea,meanPostBOLDoverRep(:,i),'color',[.5 .5 .5],...
            'linesmooth','on','linewidth',2)    
        title({'Post',['BOLD=\Sigma FR of ',num2str(numNeuroInaVoxPOST),... 
            'sampled neurons for direction ',num2str(unqDir(i)),' degrees']},...
            'fontsize',14)
        xlabel('Voxels by preferred directions (degrees)','fontsize',14)
        ylabel('Simulated Bold response','fontsize',14)
        xlim([0 numVoxInPosteriorArea])
        ylim([0.99*min(meanPostBOLDoverRep(:)) 1.01*max(meanPostBOLDoverRep(:))])
        if strcmp(sortVoxels,'sortVoxels=on')==1
            set(gca,'fontsize',14,'xtick',1:10:numVoxInPosteriorArea,'xticklabel',...
                VoxPrefDirPost(1:10:end))
        end
        box off
        
        %Prior
        subplot(313)
        plot(1:1:numVoxsPrior,meanPriorBOLDoverRep(:,i),'color',[.5 .5 .5],...
            'linesmooth','on','linewidth',2)   
        title({'Prior',['BOLD=\Sigma FR of ',num2str(numNeuroInVoxPRIOR),...
            'sampled neurons for direction ',num2str(unqDir(i)),' degrees']},...
            'fontsize',14)
        xlabel('Voxels by preferred directions (degrees)','fontsize',14)
        ylabel('Simulated Bold response','fontsize',14)
        xlim([0 numVoxsPrior])
        ylim([1.01*min(meanPriorBOLDoverRep(:)) 0.99*max(meanPriorBOLDoverRep(:))])
        if strcmp(sortVoxels,'sortVoxels=on')==1
            set(gca,'fontsize',14,'xtick',1:10:numVoxsPrior,'xticklabel',...
                VoxRespAmpPrior(1:10:end))
        end
        box off        
        drawnow
    end
end

%Distribution of voxels preferred directions for llh and Post. We expect a
%uniform representation of preferred directions across voxels.
if strcmp(sortVoxels,'sortVoxels=on')==1
    if strcmp(OptiondisplayBOLD,'displayBOLD=on')==1
        
        %print
        fprintf('%s','Calculating distribution of voxels preferred directions...')
        
        %draw
        figure('color','w')
        set(gcf,'position',...
            [0.4*screen(3) 0.05*screen(4) 0.4*[screen(3) screen(4)]],...
            'color','w')
        %llh
        subplot(121)
        hist(VoxPrefDirllh)
        set(get(gca,'child'),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
        box off
        title('llh voxels')
        xlabel('Voxels preferred directions (degrees)')
        ylabel('Number of voxels')
        box off

        %posterior
        subplot(122)
        hist(VoxPrefDirPost)
        set(get(gca,'child'),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
        title('Post voxels')
        xlabel('Voxels preferred directions (degrees)')
        ylabel('Number of voxels')
        box off
        fprintf('%s \n','done')
    end
end

%the End
fprintf('%s \n','(NNwithXGain2) Bayesian network...done') 




%DISPLAYS
%--------

%display sensory neurons responses
function display3(OptiondisplayNetwork,output,numTuningInLLH,i)

if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    figure(1)
    set(gcf,'position',...
        [0*screen([3 4]) 0.4*screen(3) screen(4)],...
        'color','w')
    clf
    subplot(531);
    set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
    hold all
    plot(output.sensNeurTuning,'linesmoothing','on','linewidth',3)
    
    %indicate currently displayed motion direction
    motionspace = 1:1:360;
    pos = motionspace(output.readoutTunings==output.motionDir(i));
    
    %indicate each neuron's activity at this direction
    plot(pos(ones(1,numTuningInLLH)),output.sensNeurFIRING(:,i),'.','markersize',20)
    xlim([0 360])
    xlabel('Motion direction (degrees)','fontsize',14)
    ylabel('Firing rate (sp/sec)','fontsize',14)
    title('llh direction-tuned neurons','fontsize',14,'fontweight','bold')
    box off
end
    
%display prior weights
function display4(OptiondisplayNetwork,output)

if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    subplot(539);
    plot(output.weightPRIOR,'color',[1 0 0],'linesmoothing','on',...
        'linewidth',3)
    xlim([0 360])
    xlabel('Motion direction (degrees)','fontsize',14)
    ylabel('Weight','fontsize',14)
    title({'Prior weight (hypothetical)','logprior(\theta_k)'},...
        'fontsize',14,'fontweight','bold','color','r')
    box off
end
    
%display sensory tuning weights
function display5(OptiondisplayNetwork,output)

if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    subplot(534)
    set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
    plot(output.readoutTunings,output.weightSensPreSynTuning,'linesmoothing','on','linewidth',3)
    xlim([0 360])
    title({'Connections weights','(logf_j)'},...
        'fontsize',14,'fontweight','bold')
    xlabel('Motion direction (degrees)','fontsize',14)
    ylabel('Weight','fontsize',14)
    box off
end

%display sensory tuning weights
function display6(OptiondisplayNetwork,output,i)
if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    subplot(537);
    plot(output.readoutTunings,output.weightLLH(:,i),'k','linesmoothing','on','linewidth',3)
    xlim([0 360])
    ylim([1.02*min(output.weightLLH(:,i)) 0.98*max(output.weightLLH(:,i))])
    title({'Pooled presynaptic inputs','(logllh(\theta_k)=\Sigma_j r_j*logf_j(\theta_k))'},...
        'fontsize',14,'fontweight','bold')
    xlabel({'Readout neurons','by readout direction (degrees)'},...
        'fontsize',14)
    ylabel('Weight','fontsize',14)
    box off
end

%display readout neurons responses
function display7(OptiondisplayNetwork,output,i)

if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    subplot(5,3,11);
    plot(output.readoutTunings,output.rPostk(:,i),'k','linesmoothing','on','linewidth',3)
    xlim([0 360])
    ylim([0.75*min(output.rPostk(:,i)) 1.5*max(output.rPostk(:,i))])
    title({'Post neurons responses','logllh(\theta_k)+logprior'},...
        'fontsize',14,'fontweight','bold')
    xlabel({'Readout neurons','by readout direction (degrees)'},...
        'fontsize',14)
    ylabel('Firing rate (sp/sec)','fontsize',14)
    box off
end


%Visualize prior and llh weights at the same scale
function display8(OptiondisplayNetwork,output,i)

if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    
    sp = subplot(538);
    hold all
    plot(output.readoutTunings,output.weightPRIOR,'color',[1 0 0],'linesmoothing','on',...
        'linewidth',3)
    plot(output.readoutTunings,output.weightLLH(:,i),'color','k','linesmoothing','on','linewidth',3)
    xlim([0 360])
    xlabel('Motion direction (degrees)','fontsize',14)
    ylabel('Weight','fontsize',14)
    title({'logPrior and logllh at same scale'},...
        'fontsize',14,'fontweight','bold','color','k')
    box off
    
    %scaling
    l = output.weightPRIOR;
    minY = min([l(:);output.weightLLH(:)]);
    maxY = max([l(:);output.weightLLH(:)]);
    set(sp,'ylim',[minY-0.1*(maxY-minY) maxY+0.1*(maxY-minY)])
end

%readout motion direction
function display9(OptiondisplayNetwork,output,i)

if strcmp(OptiondisplayNetwork,'displayNetwork=on')==1
    subplot(5,3,14)
    drawVectors(output.motionDir(i),1)
    hold on
    drawVectors(output.output.estimate(i),1)
    hold on
    drawVectors(output.PriorMean,1)
    coor = polar2cartesian(output.estimate(i),2.1);
    
    %text
    text(coor(1),coor(2),num2str(output.estimate(i)),'color','r','fontsize',20,...
        'fontweight','bold');
end
    