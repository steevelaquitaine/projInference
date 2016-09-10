
      %Author: Steeve Laquitaine
        %date: 140421
   %file name: getInstancesPerScan.m
       %usage:getInstancesPerScan
     %purpose: Get instances of BOLD for each scan of the fMRI session

%description: data are preprocessed: detrending + high pass filtering via
%concatenation function applied to individual scans
     
function [i,iscored,stimVal,pSettings,myROIname]=getInstancesPerScan(myROIname,...
    myVariable,...
    curGroup,...
    scanNumFromDescription,...
    taskNum,...
    phaseNum,...
    startLag,...
    blockLen,...
    pSettings)

%check that startLag is an integer
if (startLag~=floor(startLag))
    fprintf('\n %12s \n','(getInstancesPerScan) startLag must be an integer')
    keyboard
end

%Set parameters of the data we want to collect
%ROI: e.g., Matrix of pooled left and right MT ('allMT',140410) contains
%67 voxels (col) and 288 trials (8 scans * 36 trials (=6 dir* 3rep* 2coh)).
%2 trials were discarded at the start of each scan.
%Variable: To get variable names use "getTaskVarnames(v)" after calling
%viewSet.
%references
%http://gru.brain.riken.jp/doku.php/mrtools/mrtoolssinglepage?s[]=concattseries
v=newView;

%get the data for each scan
v=viewSet(v,'curGroup',curGroup);
scanNum=viewGet(v,'scanNumFromDescription',scanNumFromDescription,v.curGroup);
for j=1:numel(scanNum)
    v=newView;
    v=viewSet(v,'curGroup',curGroup);
    v=viewSet(v,'curScan',scanNum(j));
    
    %detrend and high-pass filter the scans you select, get rid of any
    %junk frames, do a de-trending step using a hi-pass filter if you
    %specify that, and then concatenate the resulting time series into a
    %single time series that will be saved into a new group called
    %Concatenation. The mat files are saved in folder "concatenation".
    [v,stim]=concatTSeries(v,[],'justGetParams=1','defaultParams=1',...
        'scanList',scanNum(j));
    concatTSeries(v,stim);
    fprintf('\n %12s \n',['Scan ',num2str(scanNum(j)),'(getInstancesPerScan) Concatenation...done'])
    
    %update group and scan to concatenated data
    v=newView;
    v=viewSet(v,'curGroup',stim.newGroupName);
    nScans=viewGet(v,'nScans');
    v=viewSet(v,'curScan',nScans);
    
    %extract data
    myROI=loadROITSeries(v,myROIname);
    myROI=getInstanceMatrix(v,myVariable,myROI,taskNum,phaseNum,'n=inf',startLag,blockLen);
    myROI=myROI{1};
    i{j}=myROI.classify.instances;
    stimValsUnq(:,j)=unique(myROI.classify.paramVals);
    stimVal{j}=myROI.classify.paramVals;
end
fprintf('\n %12s \n','(getInstancesPerScan) Get BOLD instances...done')



%5. Z-score i.
%Matrix must be one row for each instance, one column for each voxel
%so d=nxk where n is number of instances and m is number of voxels
%(instances (i.e., voxels responses at each trial) are z-scored across 
%trials for each voxel independently of the classes (i.e.,coherence or 
%motion direction).i.e., mean and std voxels responses ?calculated for 
%each voxel (vector of m means) and substracted (/divided) to each 
%instance of voxel response.
%If z-scoring works fine std and mean(iMatrix)=1 and 0 respectively.
[iscored,~,pSettings]=preprocessInstancesAcrossScan(i,pSettings);
fprintf('\n %12s \n','(getInstancesPerScan) Z-scoring of BOLD instances...done')

