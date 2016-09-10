
%runExperiment
%
%     author: steeve laquitaine
%       date: 140430
%    purpose: Wrapper file to run scan experiment
%
%Description: 10 scans

%first prior
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan01_Prior0_74848')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan02_Prior0_74848')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan03_Prior0_74848')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan04_Prior0_74848')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan05_Prior0_74848')

%second prior
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan06_Prior2_7714')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan07_Prior2_7714')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan08_Prior2_7714')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan09_Prior2_7714')
SLtaskLocEst('fMRI','displayName=testPsyphy','params=taskScan10_Prior2_7714')
