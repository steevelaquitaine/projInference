

%slgetDataBlockSequence.m
%
%
%
%author: steeve laquitaine
%date: 160304
%purpose:  get estimates and true directions for 100st and last trials of 
%          a previous and next prior block.
%          e.g., previous prior is 80 deg and next is 40 deg. The function
%          gets the Nst and Nlast trials of 80 deg block and the same for
%          40 deg prior (for displayed motion direction and estimates).
%
% usage:
%
%       o = slgetDataBlockSequence(d,80,10,100)
%
%Inputs:
%
%       d: databank from "SLMakeDatabank"
%   prevP: prior std previous block
%   nextP: prior std next block
%       N: number of first and last trials of each block
%
%outputs:
%
%     o.dispPrevBNst : 
%     o.estPrevBNst  :
%     o.dispPrevBNend:
%     o.estPrevBNend :
%     o.dispNextBNst :
%     o.estNextBNst  :
%     o.dispNextBNend:
%     o.estNextBNend :
%
%Description: All trials are used for the linear fit (across subjects).
%             Estimates are never averaged over any condition.


function o = slgetDataBlockSequence(d,prevP,nextP,N)

tic

%get Nst and last trials each block of the sequence
ixThisPrior = find(d.Pstd==prevP);
%get prev and next block id with
%this sequence. You can check visually:
%plot(d.Pstd); hline(80,'r'); hline(40,'r')

%when last trial of the block is the last trial of the database
%stop
if ixThisPrior(end) == length(d.Pstd)
    ixThisPrior(end) = [];
end

%When previous and next priors are the same
if prevP == nextP
    %get run ids
    thisPRunId = d.run_id(ixThisPrior);
    %find contiguous runs
    ct = 0;
    for i = 1 : length(thisPRunId)-1        
        if thisPRunId(i+1)==thisPRunId(i)+1
            ct = ct+1;
            %previous block id
            prevB_id(ct) = thisPRunId(i);
            %next block id            
            nextB_id(ct) = thisPRunId(i)+1;
        end
    end
else    
    %next blocks starting index
    a = ixThisPrior(d.Pstd(ixThisPrior+1)==nextP)+1;
    nextB_id = d.run_id(a);
    prevB_id = d.run_id(a-1);
end

%get all trials prev and next blocks
numSeq = length(nextB_id);
dispPrevBNst = [];
dispPrevBNend = [];
estPrevBNst = [];
estPrevBNend = [];
dispNextBNst = [];
dispNextBNend = [];
estNextBNst = [];
estNextBNend = [];

%extract data
%run_id = d.run_id;
run_id = d.run_id;
stimFeatureDeg = d.stimFeatureDeg;
estimatesDeg = d.estimatesDeg;

%get previous and next block trial indices in sequence
tic 
for i = 1 : numSeq
    prevthisSeq(:,i) = run_id == prevB_id(i);
    nxtthisSeq(:,i) = run_id == nextB_id(i);  
end

for i = 1 : numSeq    
    %get displayed and estimates prev block
    %     prevthisSeq = run_id == prevB_id(i);
    dispPrevB{i} = stimFeatureDeg(prevthisSeq(:,i));
    estPrevB{i} = estimatesDeg(prevthisSeq(:,i));
    
    %get displayed and estimates next block
    dispNextB{i} = stimFeatureDeg(nxtthisSeq(:,i));
    estNextB{i} = estimatesDeg(nxtthisSeq(:,i));
    
    %%--re-express all points to its angle that has the closest distance to the
    %%true direction in the space linear. This change the position of only 
    %a few edge points (3 or 4) and allows to estimate the slope between
    %estimates and displayed directions
    estPrevB{i} = slGetAngleYAsclosestToX(dispPrevB{i},estPrevB{i});    
    estNextB{i} = slGetAngleYAsclosestToX(dispNextB{i},estNextB{i});    

    %%--get first and last N trials for each block
    %prev. Fill in each variable with the first and last N trials
    %of each block
    dispPrevBNst = [dispPrevBNst;dispPrevB{i}(1:N)]; 
    dispPrevBNend = [dispPrevBNend; dispPrevB{i}(end-N+1:end)];
    estPrevBNst = [estPrevBNst; estPrevB{i}(1:N)]; 
    estPrevBNend = [estPrevBNend; estPrevB{i}(end-N+1:end)];
    %next
    dispNextBNst = [dispNextBNst; dispNextB{i}(1:N)]; 
    dispNextBNend = [dispNextBNend; dispNextB{i}(end-N+1:end)];
    estNextBNst = [estNextBNst; estNextB{i}(1:N)]; 
    estNextBNend = [estNextBNend; estNextB{i}(end-N+1:end)];
end


%%--get first and last N trials for each block
%prev. Fill in each variable with the first and last N trials
%of each block
o.dispPrevBNst = dispPrevBNst;
o.dispPrevBNend = dispPrevBNend;
o.estPrevBNst = estPrevBNst;
o.estPrevBNend = estPrevBNend;
o.dispNextBNst = dispNextBNst;
o.dispNextBNend = dispNextBNend;
o.estNextBNst = estNextBNst;
o.estNextBNend = estNextBNend;

elapsed = toc;
display(['Took ' num2str(elapsed) ' sec'])



