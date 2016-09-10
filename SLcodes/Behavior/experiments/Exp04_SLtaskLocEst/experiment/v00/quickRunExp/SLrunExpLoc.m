
%SLrunExpLoc
%
%     author: steeve laquitaine
%       date: 140430
%    purpose: Wrapper file to run scan experiment
%
%Description: 5 sessions of 5 run of psychophysics training then 5 fMRI exp. of 10 scans.



%Training
%--------

%Sess 1
dispName='displayName=VPixx';
%dispName = 'displayName=testPsyphy';

SLtaskLocEstTestGamma('psychophysics',dispName,'params=training_Prior0_74848')%80 deg
SLtaskLocEstTestGamma('psychophysics',dispName,'params=training_Prior2_7714')%40 deg
SLtaskLocEstTestGamma('psychophysics',dispName,'params=training_Prior0_74848')%80 deg
SLtaskLocEstTestGamma('psychophysics',dispName,'params=training_Prior2_7714')%40 deg





%fMRI
%----

% dispName='displayName=VPixx';
dispName = 'displayName=testPsyphy';
% dispName='displayName=fMRIproj32';

%session 1
%80 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan05_Prior0_74848')

%40 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan06_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan07_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan08_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan09_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan10_Prior2_7714')


%session 2
%40 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan06_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan07_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan08_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan09_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan10_Prior2_7714')

%80 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan05_Prior0_74848')


%3
%80 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan05_Prior0_74848')

%40 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan06_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan07_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan08_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan09_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan10_Prior2_7714')


%4
%40 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan06_Prior2_7714').

%80 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan05_Prior0_74848')


%4
%80 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan05_Prior0_74848')

%40 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan06_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan07_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan08_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan09_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan10_Prior2_7714')


%5
%40 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan06_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan07_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan08_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan09_Prior2_7714')
SLtaskLocEst('fMRI',dispName,'params=taskScan10_Prior2_7714')

%80 deg prior
SLtaskLocEst('fMRI',dispName,'params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI',dispName,'params=taskScan05_Prior0_74848')





