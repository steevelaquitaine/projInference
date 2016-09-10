
     %Author: Steeve Laquitaine
       %date: 140416 - last modification: 140421
      %usage: 
        %filename=matlab.desktop.editor.getActiveFilename
        %SLBackup(filename)
    %purpose: backup matlab .m file
%I would like to save the results. As my code is constantly changing 
%I would like to save .m Files together with the results. Just as
%a backup if I have to see exactly what I was doing.


function [mfilename,comfilename,matfilename]=SLBackup(filename)

%get the active .m file
source=filename;

%get its parent folder
%ParentFolder=fileparts(source);

%create backup directory
%backupdir=[ParentFolder,'/logs'];

%Backup file in a .txt file named with the date and time in the parent folder
%get integers only
theClock=fix(clock);

%make time readable
writableClock=[];
for j=1:length(clock)
    writableClock=[num2str(writableClock),'_',num2str(theClock(j))];
end

%append log to it to now its a backup
writableClockmFile=['_mfileLog',writableClock];

%write name new (backup) file
mfilename=[filename(1:end-2),writableClockmFile,'.txt'];

%backup
copyfile(source,mfilename)
fprintf('\n %s \n','.m file has been backed up...')

%backup commands
writableClockCom=['_CommandLog',writableClock];
comfilename=[filename(1:end-2),writableClockCom,'.txt'];
diary(comfilename)

%backup filename for .mat file
writableClockMat=['_data',writableClock];
matfilename=[filename(1:end-2),writableClockMat,'.mat'];


