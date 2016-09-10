
% author: steeve laquitaine
%   date: 150504
%purpose: initialize the parameters for each data file and get them ready
%         for analysis with "analyses" code.
%
%usage:
%
%         SLinitAnalyses


function SLinitAnalyses

%set subject
subject = input('Set subject (e.g., sub01) : ','s');

%choose analysis directory where data files will be moved
destinationDir = uigetdir(cd,'Pick the data directory where data will be moved (e.g.,/Users/steeve/data/dataPsychophy/Exp04_myProject)');
destinationDir = destinationDir;

%move to destination directory
cd(destinationDir)
mkdir 'data'
mkdir(subject);
cd(subject)

%choose data files
[datafiles,datadir] = uigetfile({'*.mat'},'Select data files','MultiSelect','on');

%move files
if iscell(datafiles) %many files
    for i = 1: length(datafiles)
        movefile([datadir datafiles{i}],[destinationDir,'/',subject])
    end
else
    movefile([datadir datafiles],[destinationDir,'/',subject])
end
fprintf('------------------------------------------- \n')
fprintf('(SLinitAnalyses) Your files have been moved :  \n')

%print files name
if iscell(datafiles) %many files
    for i = 1 : length(datafiles)
        fprintf([datafiles{i},'\n'])
    end
else
    fprintf([datafiles,'\n'])
end
fprintf(['to : ',destinationDir],' \n') %destination
fprintf('------------------------------------------- \n')

%copy and rename files 

fprintf('------------------------------------------- \n')
fprintf('Please set task info for each file ...\n')
fprintf('------------------------------------------- \n')
if iscell(datafiles) %many files
    for i = 1 : length(datafiles)
        
        %display file name
        fprintf(['Info for file : ', datafiles{i},'\n'])
        
        %set info task conditions
        exp = input('Set the experiment (e.g., 04) : ','s');
        sess = input('Set the session (e.g., 01) : ','s');
        run = input('Set the run (e.g., 01) : ','s');
        priorStrngth = input('Set the prior std (e.g., 010) : ','s');
        priorMean = input('Set the prior mean (e.g., 225) : ','s');
        stimStrngth = input('Set the stimulus strength (e.g., coherence : coh006_012 or contrast : con001_1) : ','s');
        nSampleStim = input('Set the number of sample stimuli (e.g., dir36 or loc6: ','s');
        Thedate = input('Set the date (e.g., 130504) : ','s');
        
        datafname = ['steeve_exp',exp,'_data_',subject,'_sess',sess,'_run',run,'_Pstd',priorStrngth,'_mean',priorMean,'_',stimStrngth,'_',nSampleStim,'_',Thedate,'.mat'];
        
        copyfile(datafiles{i},datafname)
        
        %display file name
        fprintf(['The file ', datafname,' has been created \n \n'])
        
        
    end
else 
    %display file name
        fprintf(['Info for file : ', datafiles,'\n'])
        
        %set info task conditions
        exp = input('Set the experiment (e.g., 04) : ','s');
        sess = input('Set the session (e.g., 01) : ','s');
        run = input('Set the run (e.g., 01) : ','s');
        priorStrngth = input('Set the prior std (e.g., 010) : ','s');
        priorMean = input('Set the prior mean (e.g., 225) : ','s');
        stimStrngth = input('Set the stimulus strength (e.g., coherence : coh006_012 or contrast : con001_1) : ','s');
        nSampleStim = input('Set the number of sample stimuli (e.g., dir36 or loc6: ','s');
        Thedate = input('Set the date (e.g., 130504) : ','s');
        datafname = ['steeve_exp',exp,'_data_',subject,'_sess',sess,'_run',run,'_Pstd',priorStrngth,'_mean',priorMean,'_',stimStrngth,'_',nSampleStim,'_',Thedate,'.mat'];
        copyfile(datafiles,datafname)
        
        %display file name
        fprintf(['The file ', datafname,' has been created \n \n'])
end
