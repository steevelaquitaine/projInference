
% Steeve laquitaine 19032009

% need log files (with their actual name)in the dir
% rename each file as : f1.......fn
% give for each log files: 
     % R: local action


clear all
r = dir;
r(1:2,:) =[];

for i = 1:size(r); % select the files
        %str = r(i).name;
        Name{i}=r(i).name;
        %date{i}= str(1:8);
        %proba{i}=str(10:13);
        %session{i}=str(15:18);
        
        fid = fopen(r(i).name,'r'); % id the files
        f = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',[18 inf]);
        f=f';
        eval(['f' int2str(i) '=f';]);
        
        %success trials
        t=find(f(:,2)==1); %trial number
        
        %action
        R=f(t,3); % success trials
        eval(['R' int2str(i) '=R';]);
        clear t

 
        R2=R; 
        R2(length(R2)+1:900)=NaN;
        reward(:,i)=R2;
        clear('R','R2')
        
end  


clearex('reward','Name')


