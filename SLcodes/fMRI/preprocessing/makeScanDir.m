
% author: Steeve Laquitaine
%purpose: get tar file from NIMS website, process the content files and
%         sort them into Anatomy and Raw/Tseries folders for mrLoadRet
%
%cd to your .tar file's directory

function makeScanDir(tarfilename)

%decompress file
!tar -xvf tarfilename

%move last folder to s0301150324

%uncompress files with each branch folders
%-rename files with branch folders nam
%-create Anatomy folder
%-move anatomy files to Anatomy folder
%-delete branch folder (now empty)
%