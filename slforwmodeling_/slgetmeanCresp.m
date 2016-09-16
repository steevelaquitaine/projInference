
%slgetmeanCresp.m 
%
%
% author: steeve laquitaine
%purpose: reposition channel responses as a function of the distance of
%         their preference to the displayed direction which is always 
%         position at the center
%
%
%inputs : 
%
%       Cresp :
%        svec :
%       phi_k : is a Nv directions x 1 column vector
%      center :
%
%
%output
%
%       fm.meantestCresp : 1 x Nk channels vector of average channel responses 
%                          where channel responses have been repositionned
%                          relative to the distance of their preference to the
%                          displayed direction (always set at the "center" 
%                          argument position)
%
%         fm.alltestCresp : Ni instances x Nk channels matrix of all channel 
%                          responses re-positioned
%
%
%reference : inspired from Justin Gardner, cinvor code

function fm = slgetmeanCresp(Cresp,svec,phi_k,center)

%check that channel preferences match displayed directions
%phi_k is a Nv directions x 1 column vector
phi_k = SLmakeColumn(phi_k);
if ~isequal(unique(svec),phi_k)
  fprintf('%s \n','(slgetmeanCresp) Channel preferences do not match displayed directions: centering channel responses with preferences matching the displayed direction is not possible');
end

%find center position
centerIdx = find(phi_k==center);

%stop if no channel preference match the center
if isempty(centerIdx)
  fprintf('%s %i \n','(slgetmeanCresp) One channel should be centered on', center);
  keyboard;
end

%loop over channels by direction preference
fm.alltestCresp  = [];
for i = 1 : length(phi_k)
  %get all channels' responses for the trials where the displayed
  %direction matched this channel preference
  theseRespIdx = svec==phi_k(i);
  theseStimResp = Cresp(theseRespIdx,:);  
  %shift the column position of this channel to the center
  %by how far it's preference is from the displayed direction
  fm.alltestCresp = [fm.alltestCresp ; circshift(theseStimResp,[0 centerIdx-i])];
end
%average the displayed-direction-re-centered channel responses
fm.meantestCresp = mean(fm.alltestCresp ,1);
