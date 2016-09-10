
% steeve laquitaine 23082009

% Transform Bachoice and GoodChoice  ts from nex into scalar action (1 and 0) 

Badchoice(:,2)=0; 
GoodChoice(:,2)=1 

choice=[Badchoice;GoodChoice]; 
[d1,d2] = sort(choice(:,1));

% choices
action=choice(d2,:)

%clearex('action','Badchoice','GoodChoice');