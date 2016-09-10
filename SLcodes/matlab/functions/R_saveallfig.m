
clear all
r = dir;
r(1:2,:) =[];

for i = 1:size(r); % select the files
        %str = r(i).name;
        %Name{i}=r(i).name;
        %date{i}= str(1:8);
        %proba{i}=str(10:13);
        %session{i}=str(15:18);
        
        fid = fopen(r(i).name,'r'); % id the files
        f = fscanf(fid,'%g',[150 inf]);
        f=f';
        
        y=surf(f)
        shading interp 
        H=get(0,'children'); %get(root handle,'property');handles of children of the root object (one figure)
        saveas(H,['shyashin' num2str(i) '.emf']) % file can be '.fig' or else
        clear H
end

 