
% steeve laquitaine 23082009

% Transform Bachoice and GoodChoice ts from nex into scalar action (1 and 0) 


noreward(:,2)=0
reward(:,2)=1

rnfrcmt=[noreward;reward]
[d1,d2] = sort(rnfrcmt(:,1));

% reinforcement
reinf=rnfrcmt(d2,:)
reinf=reinf(:,2);

%clearex('reinf','noreward','reward');