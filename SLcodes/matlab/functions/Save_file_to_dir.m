



% SAVE IN DIRECTORY

 file=whos('PCTH*');
 
    for i=1:size(file,1);
        file_name{i}=file(i).name;
        save (file_name{i},file(i).name);
    end