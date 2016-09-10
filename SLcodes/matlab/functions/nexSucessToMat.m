
% steeve laquitaine 20090914

% need Int_Success_Trials & trial_begin
% compute Successful trials in a vector where 1 is for successful trials 0
% for motor error trials 


Trials=[trial_begin zeros(size(trial_begin,1),1)]
for i =1: size(Int_Success_Trials,1)
        [c(i)]=find(trial_begin==Int_Success_Trials(i,1)) 
end
Trials(c,2)=1;
S=Trials(:,2);