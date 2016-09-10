

%slfmriStckAddSwitchingVar.m
%
% author: steeve laquitaine
%purpose: calculate and add switching variable to the stim files of each
%         fmri session of the project


%set fMRI sessions
o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';
o.sessPath{1} = [o.rootpath '/s02520150814'];
o.sessPath{2} = [o.rootpath '/s002520150923'];
o.sessPath{3} = [o.rootpath '/s002520150925'];

%set group
o.myGroup = 'Concatenation';

%set task
o.taskNum = 2;

%calculate and add switching variable to each session
for sess = 1 : length(o.sessPath)
    fprintf('%s%s \n','SESSION:',o.sessPath{1})
    fprintf('%s \n','--------------------------------------------')    
    o.myPath = o.sessPath{sess};
    tmp = slmakeSwitchingVar(o);
    o.behBySess{sess} = tmp;
end
