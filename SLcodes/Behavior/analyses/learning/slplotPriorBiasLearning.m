


%%get data
% databank = SLMakedatabank({'sub01','sub02'},'dataPath',...
%     '~/data/dataPsychophy/proj01_priorStrength/',...
%     'experiment','vonMisesPrior');
% 
% slplotPriorBiasLearning(databank)


function slplotPriorBiasLearning(databank)

d = databank;

%get variables
feat = d.stimFeatureDeg;

%get data
estimates = d.estimatesDeg;



