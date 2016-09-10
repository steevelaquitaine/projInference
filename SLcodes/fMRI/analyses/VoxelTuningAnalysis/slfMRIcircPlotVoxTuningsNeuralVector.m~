
%slfMRIcircPlotVoxTuningsNeuralVector.m



function slfMRIcircPlotVoxTuningsNeuralVector

cd('/Volumes/DroboBKUP/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/myRandomCoh/accAtTimeleaveOneOutfisherbalancByRemovI/outsideBrain02/')
load voxRawTunings.mat 
%%
%loop over voxels
for i = 1:53

    figure('color','w')
    plot(st1bkp{i}.conditions(:,1),st1bkp{i}.mean)
    a = SLcircWeightedMeanStd(st1bkp{i}.conditions(:,1)',st1bkp{i}.mean/nansum(st1bkp{i}.mean)); a.deg.mean
    SLcircHist([15 85 155 295],st1bkp{i}.mean,[1 0 0],40)
    hold on; drawVectors(a.deg.mean,200)

end

