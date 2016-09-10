
%slSimBayesLearning.m
%
%
%     author: steeve laquitaine
%       date: 150919
%    purpose: Learning prior expectation about stimulus distribution with Bayesian filtering
%
%Description: The past posterior is the prior for the current trial. It summarizes 
%			  the current trial prior belief built on past sequence of evidence:
%
%			  P(M|D) ~=  P(D|M) * P(M)
%			  P(M|D1,D2...) ~=  P(D1,D2|M) * P(M) = P(D2|M) * P(D1|M) * P(M)
%
%			  where P(M) is our initial prior
%	   		  P(D1|M) * P(M) the updated prior after D1 (a posterior)
%
%
%  reference:  p9, Bayesian Brain by Kenji Doya et al.

kappa             = 10; 								     %sensory noise
kappaS            = 80;								 		 %prior strength
stimulusSpace     = 1:1:360;
prior             = vmPdfs(stimulusSpace,225,kappaS,'norm'); %stimulus statistics (correct prior)
stimuli           = randsample(stimulusSpace,100,'true',prior); %stimuli time series
evidenceSpace     = 1:1:360;
evidenceDensities = vmPdfs(evidenceSpace,stimuli,40,'norm'); %evidence density given stimulus

%learn prior based on evidence
priorSubj = ones(360,1)*1/360; %initial subjective prior

for i = 1 : length(stimuli)
	
	%observed evidence from stim.
	evidence(i)      = randsample(evidenceSpace,1,'true',evidenceDensities(:,i)); %evidence from stim

	%infer back stimulus from evidence
	like(:,i)        = vmPdfs(stimulusSpace,evidence(i),kappa,'norm') ;      	  
	posterior(:,i)   = like(:,i).*priorSubj(:,i)/sum(like(:,i).*priorSubj(:,i));

	%update prior about stim stats
	priorSubj(:,i+1) = posterior(:,i); 

end 

imagesc(priorSubj)


%% long tailed likelihood
%long-tailed likelihood produces long-tailed prior.

kappa             = 10; 								     %sensory noise
kappaS            = 80;								 		 %prior strength
stimulusSpace     = 1:1:360;
prior             = vmPdfs(stimulusSpace,225,kappaS,'norm'); %stimulus statistics (correct prior)
stimuli           = randsample(stimulusSpace,100,'true',prior); %stimuli time series
evidenceSpace     = 1:1:360;
evidenceDensities = vmPdfs(evidenceSpace,stimuli,40,'norm'); %evidence density given stimulus

%learn prior based on evidence
priorSubj = ones(360,1)*1/360; %initial subjective prior
t = 0.3;

for i = 1 : length(stimuli)
	
	%observed evidence from stim.
	evidence(i)      = randsample(evidenceSpace,1,'true',evidenceDensities(:,i)); %evidence from stim

	%infer back stimulus from evidence
	like(:,i)        = (1-t)*vmPdfs(stimulusSpace,evidence(i),kappa,'norm') + t;      	  
	posterior(:,i)   = like(:,i).*priorSubj(:,i)/sum(like(:,i).*priorSubj(:,i));

	%update prior about stim stats
	priorSubj(:,i+1) = posterior(:,i); 

	plot(priorSubj(:,i))
	drawnow
	pause(1)

end 

%imagesc(priorSubj)
