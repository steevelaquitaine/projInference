
% Steeve laquitaine 19032009

% need log files (with their actual name)in the dir
% rename each file as : f1.......fn
% give for each log files: 
    % p_a : smoothed probability of A1 for window of bin trials
    % A: local action


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
        f = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',[21 inf]);
        f=f';
        eval(['f' int2str(i) '=f';]);
        
        %success trials
        t=find(f(:,2)==1); %trial number
        
        %action
        A=f(t,4); % success trials
        eval(['A' int2str(i) '=A';]);
        clear t

 
        A2=A; 
        A2(length(A2)+1:900)=NaN;
        action(:,i)=A2;
        clear('A','A2')
        
end  


clearex('action','Name')


