
% Steeve laquitaine 14062009


clear all
r = dir;
r(1:2,:) =[];

for i = 1:size(r); % select the files
        str = r(i).name;
        Name{i}=r(i).name;
        date{i}= str(1:8);
        proba{i}=str(10:13);
        session{i}=str(15:18);
        
        fid = fopen(r(i).name,'r'); % id the files
        f = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',[18 inf]);
        f=f';
        eval(['f' int2str(i) '=f';]);
        
        % performance 
        t=find(f(:,2)==1); %trial number
        Perfo(i)=mean(f(t,4)); % success trials
        % motor success
        m_Success(i)=mean(f(:,2));     
end
        
        
        