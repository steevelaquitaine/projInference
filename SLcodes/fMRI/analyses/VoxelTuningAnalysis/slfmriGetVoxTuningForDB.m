
%slfmriGetVoxTuningForDB.m


function [st1,modes1,st2,modes2,nVox] = slfmriGetVoxTuningForDB(d,tuningType,o)

nVox = size(d.instances,2);
tic
if strcmp(tuningType,'vmfit')
    parfor voxnum = 1:nVox
        st1{voxnum} = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch=1','myRandomCoh=0.08',d,'gridsearchMean');
        modes1(voxnum) = st1{voxnum}.fitvmMode; %selectivity for sw-p
        st2{voxnum} = slfmrigetVoxTuning(voxnum,'myRandomDir','mySwitch=2','myRandomCoh=0.08',d,'gridsearchMean');
        modes2(voxnum) = st2{voxnum}.fitvmMode; %for sw-d
        vline(225,':')
    end
end

if strcmp(tuningType,'neuralvector')
    %equalize direction distributions
    %across switching variable conditions
    [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist00(d,o);
    parfor voxnum = 1:nVox        
        
        %switch to prior
        st1{voxnum} = slfmrigetSingleVoxTuningNeuralVector(voxnum,'myRandomDir','mySwitch=1',['myRandomCoh=' num2str(o.coh2plot)],d_dirBalbySw{coh_u==o.coh2plot},'gridsearchMean');
        modes1(voxnum) = st1{voxnum}.fitvmMode; 
        
        %switch to direction
        st2{voxnum} = slfmrigetSingleVoxTuningNeuralVector(voxnum,'myRandomDir','mySwitch=2',['myRandomCoh=' num2str(o.coh2plot)],d_dirBalbySw{coh_u==o.coh2plot},'gridsearchMean');
        modes2(voxnum) = st2{voxnum}.fitvmMode; 
        vline(225,':')
    end
end
toc
close all