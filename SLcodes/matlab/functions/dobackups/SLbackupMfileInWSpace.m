
%SLgetActivemFile.m
%
%authors: steeve laquitaine
%   date: 140801
%purpose: backup active mfile in the workspace in  "myMfile" variable.
%
%usage: 
%      mfilename = SLgetActivemFile
%      myMfile = SLbackupMfileInWSpace(mfilename )

function myMfile = SLbackupMfileInWSpace(mfilename)

%write in in worspace
fid = fopen(mfilename, 'r');       
myMfile = fread(fid, inf, 'uint8=>char')';
