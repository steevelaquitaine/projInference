# slCodes
________________________________________________________________________
**Author**: Steeve Laquitaine, <steeve@stanford.edu>

**Description**: Library of codes to study how humans use knowledge of motion direction statistics 
to improve motion direction estimation.

**Overview**:
- behavioral experiments
- fMRI experiments 
- and data analyses.

**Codes**:
________________________________________________________________________

##### Directories

SLexistPath


##### Database operations

[slgetdatapath](Behavior/data/slgetdatapath.m) : retrieve data path for motion direction inference project ("lab" and "home").  

[SLMakedatabank.m](Behavior/analyses/preprocessing/SLMakedatabank.m): make matrix database of motion direction estimate data, associated variables (motion coherence, direction, motion direction prior standard deviation, prior mean) and additional info such as trial number, behavioral session, run.

[SLinitFAs.m](Behavior/analyses/preprocessing/SLinitFAs.m): initialize the variables (coherence, direction, direction prior standard deviation) associated with the motion direction estimate data.  

[SLfindEmptyCells.m ](matlab/functions/dataWrangling/SLfindEmptyCells.m)  : find empty cells in a matrix

[SLfindRow.m](matlab/functions/dataWrangling/SLfindRow.m)  : find the position of a specified values of a row in a matrix.  

SLuniqpair  
SLisCellEmpty


_csv_                                 

[slcsvRead](matlab/database/slcsvRead.m) : read .csv data table  
[slcsvReadNumCol](matlab/database/slcsvReadNumCol.m) : count the number of columns in a csv file  
[slcsvReadNumRows](matlab/database/slcsvReadNumRows.m) : count the number of rows in a csv file

#### Data wrangling

[sladdZerobeforeNum2str.m](matlab/functions/dataWrangling/sladdZerobeforeNum2str.m) : add a zero to a number when converting to string '01' instead of '1';



##### Backup
SLautobackup        
SLBackup                                                
SLbackupMfileInWSpace                                   
SLsaveTightFigure

##### Text mining

SLLoadTextData                                          
SLfindword


##### Numbers
SLgetNumPrecision  
SLisinteger

##### Data transformation

[slScaleBetween0and1](matlab/functions/scaling/slScaleBetween0and1.m) : Scale data between 0 and 1

##### Geometry

SLcircle

##### Algebra

SLloglogpn  
SLdrawSigmoid  
SLdrawSoftmax

##### Linear algebra (vector operation)


SLmakeColumn  
SLmakeRow  
SLcalculateVectNorm   
SLcart2polar    
SLcircCcwVect2Angle   
SLpolar2cartesian   
SLreplicateRows   
SLeuclidDist    
drawVectors : draw vectors    
[slgetVectorStats](matlab/functions/vectors_calculus/slgetVectorStats.m) : calculate vector length and directions statistics probability distribution of vector direction, each direction average vector length and sem.  
[slsubstRowAndcolVec.m](matlab/functions/matrix/slsubstRowAndcolVec.m) : quickly substract row and col vectors x and y in the following way :
```matlab
x = 1:10; y = x';
for i = 1 : length(y)
   xminy(i,:) = x - y(i);
end
```    

[slissymetric.m](matlab/functions/matrix/slissymetric.m): Check whether matrix is symmetric.  



## Graphics
SLlinspecer   
SLBreakBar  
SLBreakPlot  
SLcircPlot  
SLcircShiftXandY4Plot      
[sldrawVecPdf](/matlab/graphics/circular/sldrawVecPdf.m) : plot vector probability distribution in a circular graph     
[sldrawCircErrorBar](/matlab/graphics/circular/sldrawCircErrorBar.m) : circular plot mean vector length and error bar by vector directions


slgetAllpolarOjects : get axis object handles for texts, lines, axes and patch    
sldeletePolarOuterLine : remove the polar thick outerline in a figure

SLConventionUp  
SLdarken  
SLplotLabelAll  
SLplotNumRowAxes  
SLdrawBar  
SLpositionFigure  
SLpositionFigures  
SLquickPlot  
SLremoveDeadSpace  
SLdrawPublishAxis  
SLerrorarea  
SLerrorbar  
SLsetFigPublishSize  
SLshadesOfGreyLinspecer  
SLgetAxes  
SLgetFigures  
SLgetLegendAxes   
SLprintProgress   
SLprintTable  
SLwhiten  

SLisXLabel: says if an axis has an xlabel (1) or not (0)                                          

slGraphWdgLines : create a wedge matrix to display with imshow or use with mglCreateTexture.m (mgl distribution)  
slGraphBlurEdge : create a blurred edge matrix to display with imshow or use with mglCreateTexture.m (mgl distribution)  

slMakeLocStimWdgAreaGrBg: location stimulus contrast is such that within wedge area region lighter than the mean luminance within wedge area see an increase in their luminance by contrast/2 and darker one a decrease by contrast/2. Gray bkg.  

slMakeLocStimWdgAreaBgFlash :location stimulus contrast is such that within wedge area region lighter than the mean luminance within wedge area see an increase in their luminance by contrast/2 and darker one a decrease by contrast/2.Flashing cloudy bg.

slMakeLocStimWdgArea :location contrast wedge stimulus is such that within wedge area region lighter than the mean luminance within wedge area see an increase in their luminance by contrast/2 and darker one a decrease by contrast/2.
 
slMakeLocStimBlurEdg: location contrast edge stimulus. Stimulus is a blurred oriented edge (hallow of light) on a gray background
 
slMakeLocStimMasked : brief contrast black and white flickering wedge displayed at different locations on gray Bg  followed briefly by lo-f or hi-f bg (black and white mask)
 
slMakeLocStimMixtCon: contrasted wedge lines location stimulus such there is a change in contrast within wedge area: increase in luminance by contrast/2 and decrease by contrast/2 relative to average bg luminance

slMakeLocStimInCloud: uniformly bright radial lines displayed at different locations superimposed on a cloudy contrasted  background


## Statistics

#### Descriptive statistics

SLmakeMean     
[slCalcMeanByRowCondEachCol](matlab/functions/proba_statistics/slCalcMeanByRowCondEachCol.m) : calculate the means over conditions labelled in the rows of an M by N matrix

SLmakeVar  

SLmakeStd   

[slMakeCI](matlab/functions/proba_statistics/slMakeCI.m) : calculate margin error, confidence interval, mean and standard deviation of a data sample (the confidence interval is given by +- gaussian criterion for your confidence (e.g. 95%) * the standard error of the mean (std/square root the sample size))

SLmakeStat: sort conditions and get linear stats (mean and std) for each  

**_SLmakeDataMeanAndStd_**  

**_slBootstrap_** : create 'nBoot' new datasets by boostrapping {x,y} data pairs and an arbitrary operation on each dataset with output stat (e.g.,mean)

**_SLbootstrapMedianStd_**
**_SLmakeDiscreteMixtureGauss_**  
**_SLmakeStat_**  
**_SLmixtureGaussian_**  
**_SLcumSumInBin_**  
**_SLPoissonPdfs_**  
**_SLGaussianPdfs_**  
**_SLgetStat_**  
**_SLzscore_** :  z-score data


####  Circular statistics
SLvectors2signedAngle  
SLmakeCircStat  
SLcircConv  
SLmakeDiscreteMixtureVM  
SLcircHist  
SLcircMeanStd  
SLcircMeanStd2  
SLcircBin  
SLCircStat  
SLcircWeightedMeanStd   
SLMixtureVM  
SLdegLin2Circ  
SLpolar  



**_SLdrawCircStat_ : draw a circular data mean and std against circular independent variables "v1" and up to two more variables "v2","v3" (e.g., dataRaw = mean motion direction estimates and v1 = true motion direction, v2 = motion coherence, v3 = direction statistical distribution strength)
SLdrawVonMisesStdofKarea  
SLfitMixtVMtoMixtGauss  
SLstatcircular   
SLKtoStdev: approximate circular distribution concentration parameter kappa in standard deviation

slGetAngleYAsclosestToX : convert angle y (deg) in the linear space so that it is expressed as the closest angle to x (in signed deg). 



#### Modeling tools

SLmakeR2 : calculate the percent of variance in a dataset that is explained by a model that is fiited to the dataset

SLmodelComparison

[slGetAIC](https://github.com/steevelaquitaine/SLcodes/blob/fb2eb6e7a8cc6949944645849a81e73240a10f90/Behavior/analyses/modeling/modelComp/slGetAIC.m) : calculate Akaike information criterion for model comparison.

SLdrawAICperCondition: plot AIC (akaike information criterion) from maximum likelihood fitting of data given a model.

slCalculateAICbyCondition : calculate AIC for a given model and selected experimental conditions.

slCalculateAICbyCondition_loop : calculate difference in AIC between a reference model (switching) and the other Bayesian models for selected experimental conditions and plot them.

slPlotModelsAICvsSwitchingAIC : calculate difference in AIC between a reference model (switching) and the other Bayesian models for all  experimental conditions together and plot them.

SLdrawModelParameters

SLgetFitParamStd

[slrecordfit](Behavior/analyses/slrecordfit.m): record fit data (parameters explored during search, objective function values, iterations) in a "history" variable.
<p align="center">
  <img src="matlab/optimization/figures/recordfit.png?raw=true" alt="Sublime's custom image"/>
</p>

[slLoadSavedFitParams](/Behavior/analyses/modeling/models/slLoadSavedFitParams.m) : load model saved fit data (best fit parameters, logl, etc ..)

#### Bayesian modeling of motion direction estimation

[SLfitBayesianModel](https://github.com/steevelaquitaine/SLcodes/blob/fb2eb6e7a8cc6949944645849a81e73240a10f90/Behavior/analyses/modeling/models/Bayesian/SLfitBayesianModel.m): Fit Bayesian models to motion direction estimation data:
- Basic Bayesian model with circular evidence distributions (Stocker et al, 2006; Girshick et al., 2011) 

```matlab
SLfitBayesianModel({'sub02'},[80 40 20 1.74 4.77 10.74 34.25 NaN 0.001 15 NaN],'experiment','vonMisesPrior','filename','datafit','MAPReadout','MaxLikelihoodFit','fminsearch');
```
- Bayesian Sampling (Carandini et al, 2011)

```matlab 
SLfitBayesianModel({'sub01'},[80 40 20 1.74 4.77 10.74 34.25 NaN 0.001 15 NaN],'experiment','vonMisesPrior','filename','datafit','SamplingReadout','MaxLikelihoodFit','fminsearch');
```

- Basic Bayesian model with cardinal priors (Girshick et al., 2011)

```matlab
%see SLfitBayesianModelExamples for example model setups
>> edit SLfitBayesianModelExamples.m 
```

SLBayesSamplingLookupTable                          
SLmakePredictionsBayesianModel  
SLdrawModelRepresentations   
SLdrawModelsPredictionCentered
SLdrawModelsPredictionHistwithDataCentered  SLgetLoglBayesianModel  
SLGirshickBayesLookupTable
slgetloglWJMfromfitp : get logl of data given parameters for the Bayesian model with long-tailed likelihood and circular prior.

slPlotPropSubSwitchingVsBayesModels :  plot proportion of subjects for which switching outperforms other models.


##### Heuristic approximation to Bayesian estimation

SLfitCompetitionModel  
SLmakePredictionsCompetitionModel


#### Mechanistic modeling of brain activity

SLmakeTuningCurves

SLdoFEandReconstruct : forward encoding modeling and reconstruction

SLdoForwardEncoding : forward encoding modeling

SLfMRIforwardModelingAnalysis

SLsortByParamSelectivity

SLsortTrainingAndTest                                  

SLgetParamSelectivity                                                                      

slfmriClassifStckS.m : classify from fMRI sessions data stacked in one big database



#### Simulations

[slsimulateBayesianModelEstimateSpace](/Behavior/analyses/modeling/simulations/slsimulateBayesianModelEstimateSpace.m) : simulate estimate space produced by Bayesian inference with circular (von Mises) vs. linear (Gaussian) likelihood and priors. Estimate space is different because in circular space evidence can be pulled from both side of the prior when evidence distribution peaks on the opposite side of the prior while in linear space from only one side. This produces smaller biases and larger variability further from prior in circular space while constant bias and variability in the linear space.

[slsimulateBayesianModelCircularEstimateSpace](/Behavior/analyses/modeling/simulations/slsimulateBayesianModelCircularEstimateSpace.m) : simulate estimate space produced by Bayesian inference with circular likelihood and priors.

[slsimulateBayesianModelLinearEstimateSpace](/Behavior/analyses/modeling/simulations/slsimulateBayesianModelLinearEstimateSpace.m) : 
simulate estimate space produced by Bayesian inference with Gaussian likelihood and priors in a linear space.

SLsimAttractor  
SLsimAveragePPCafterAdaptation  SLsimAveragePPCafterAdaptationBaseChange  
SLsimAveragePPCafterAdaptationGainChange  
SLsimAveragePPCafterAdaptationSpikeRangeChange  
SLsimBayesianPredictionsWithExpPrior  
SLsimCompetitionPredictionsWithExpPrior  
SLsimPopVoxAnalforBayesianIntegVsCompetition  
SLsimPPC  
SLsimPPCafterAdaptation  
SLsimulateAdaptedPrior  
SLsimulateBayesianModelWithData  : simulate a basic Bayesian model estimate mean and variability predictions superimposed with data from a subject in a motion direction estimation task 
slSimulateBayesianModel :  simulate a basic Bayesian model estimate mean and variability predictions in a motion direction estimation task
SLsimulateBAYESmapAndSamplingAndCompetitionPredictions  
SLsimulateBayesVsCompBimodalPrior  
SLsimulateCompetitionModel  
SLsimulateFatPTweightedPrior  
SLsimulateFatTailPrior  
SLsimulateFatTailPriorAndLlh  
SLsimulateFatTailPriorAndPoissonLLH  
SLsimulateFatTailRepresentations  
SLsimulateLaplacePrior  
SLsimulatePoissonPriorAndLLH  
SLsimulatePriorAndPoissonLLH  SLsimulatePTweightedPriorAndLLH  SLstepBystepBayesianObserver

### fMRI

[slfmriSimInstances.m](fMRI/database/slfmriSimInstances.m) :  simulate instances of Bold responses 
to a stimulus that will be use for classification.


## Machine Learning

SLfMRIclassificationAnalysis  

[slfmriClassifyByConditions](fMRI/analyses/classification/slfmriClassifyByConditions.m) : classify motion directions for different conditions (e.g., switching). Equalize the number of instances between the conditions to allow a fair comparison of the classification accuracies between conditions.
   
## Audio
SLmakeSound
                                               
                                           
## Programming
SLpauseCode  
SLgetActivemFile  
sldispPerc : display progress in a loop in percent                                                                                                                                              
 


## Experiments    
#### Set visual stimulus parameters
SLtransfDist2VisAngle  
SLtransfPixelstoCmeters  
SLvisangle2stimsize
                            
                                                                                                                                            
#### Set experiment parameters 
SLinitRunExpBimodalPrior  
SLrunExpDotDirBimGmixtPrior  
SLrunExpDotDirBimVMmixtPrior  
SLrunExpDotDirUniPrior


#### Eye tracking analysis
SLgetTaskEyeTraces :  
SLanalysesEyeMvt : plot eye positions and compare them between desired conditions with Anova.



## Behavioral analysis
______________________

##### : Motion direction estimation behavior ("proj01_priorStrength")

*Goal*: model how humans use the statistical distribution of motion direction to improve motion direction estimation 

*data preprocessing and organization workflow*  
[slworkflowRenameBehData.m](/Users/steeve_laquitaine/DropBox/myDropBox/Codes/SLcodes/Behavior/workflow/slworkflowRenameBehData.m) : preprocess the raw data from fMRI sessions and organizing them in structured path. Files are renamed with session, run, date and task conditions info to be analysed with the script "analyses.m".  

[slplotandStatsPrevDirvsPriorBias](Behavior/analyses/prevDirEffect/slplotandStatsPrevDirvsPriorBias.m) : plot the average bias toward the prior versus the previous motion direction relative to the current direction for each subject. Also tests the null hypothesis that the average overall distance between the estimate and the prior over all the trials (over conditions) is > 0 and that its lower 95% confidence is larger than 0 (bias toward prior). Also test the hypothesis that estimate is biased toward previous direction.

[slgetprevDirEffect](Behavior/analyses/prevDirEffect/slgetprevDirEffect.m) : get all the trials in which the currently displayed direction is in between the previous direction and the prior mean and plot the average estimates to test whether they are biased toward the previous direction or the prior mean.

**slwrapperPlotPriorBiasLearning** : analyse prior metalearning : the effect of previous prior block on learning in next prior block

**slAnalysisLearningScatterEarlyVsLateBias** : purpose: Scatter plot the slopes of estimates vs true motion direction for
early versus late trials pooled across runs. N Early and N late trials are sorted for each prior and pooled together Trial data are linearized as data + distances to true motion direction (now range between -180 deg and 540) then          averaged for each direction, coherence and prior condition. The means are then fitted with a line producing a slope estimate for each condition.

**slgetDataBlockSequence**  :  
[slIsPriorStrgDisNormal.m](Behavior/analyses/modeling/priorStrengths/slIsPriorStrgDisNormal.m) :check whether prior strengths are normally distributed across subjects for the Basic Bayesian observer and the switching observer p < 0.05 indicates that the null hypothesis that the data are normally distributed is rejected (a non parametric test should be used).


#### database (outside scanner training)

[slrenameFileForAnalysis](Behavior/analyses/preprocessing/slrenameFileForAnalysis.m) : rename raw stim files for routine analysis of behavioral data from training subjects before fMRI scanning.  



## fMRI analysis
________________

notes : requires cloning mgl distribution

##### database  

[loadOrReloadROITseries.m](fMRI/database/loadOrReloadROITseries.m) : load and reload ROITseries faster by saving them in a directory once loaded once.  

[slfmriGetDBoverSessions.m](fMRI/database/slfmriGetDBoverSessions.m) : retrieve database built over sessions of trial instances of voxel response to the motion stimulus with their associated variables for each trial (motion direction, coherence, subject switching to prior or evidence, timing of stimulus acquisition volumes).        

[slfmriCreateOrLoadfMRIdb.m](fMRI/database/slfmriCreateOrLoadfMRIdb.m): Create fMRI database structure "d.mat" containing Bold response instances, task variables and timing of volume acquisition from fMRI time series and stimulus "stim.mat" files or load already existing database structure "d.mat".  

[slsimfMRIdatabase](fMRI/database/slsimfMRIdatabase.m) : simulate fMRI database of trial instances voxel responses with their associated variables. Motion direction can be von Mises distributed across trials ('vmDist'); the voxels responses can be random values ('random'), or there can be 5 voxel selectivities equiprobably distributed (no arg).

[slfmriStackAndSaveSessionInstances](fMRI/database/slfmriStackAndSaveSessionInstances.m) : stack and save session instances in a structured directory


##### visualize raw database  

[slfmriVisualDatabase](fMRI/database/slfmriVisualDatabase.m) : 


##### Functional localization analyses  

[slfmriGetMotLocr2](fMRI/analyses/Localization/slfmriGetMotLocr2.m) : load the r-square of fitting brain voxels' responses with a sinewave with the same period as an on/off motion stimulus localizer. Voxels with high r-squared are considered functionnally responsive to motion.

##### full analyses

slfmriScriptERana : run fMRI event-related analysis  
SLfMRIcorrelationAnalysis :  
SLfMRIeventRelatedAnalysis :  



##### machine learning
slfmriInitClassifAnalysisTaskDotDirfMRI05 :  set parameters for fmri classification analyses  
slfmriClassifyDir_x_switch : classify motion directions from brain activity for data sorted by "switching" variable to determine whether sensory evidence is represented or not when subject switch to prior  
slfmriClassify : decoding stimulus features (e.g., motion direction) from BOLD signals  
slloadChancePermutationData : load chance permutation data for slTaskDotDirfMI05 project.  
slfmriClassify2 : decoding various features with various options from ROI voxels Bold pattern. Get Bold response
instances.



##### voxel tuning  
slfmrigetVoxTuning : plot voxel tuning over variable 0 sorted by variable 1  
[slfmriGetSingleVoxNullneuralVectorDist](fMRI/analyses/VoxelTuningAnalysis/slfmriGetSingleVoxNullneuralVectorDist.m) : calculate single vox null neural vector distribution, by permutation of responses y association to stimulus circular feature x (e.g., motion direction, orientation etc...) and voxel actual neural vector.  
[slfmriGetVoxTuningParams](fMRI/analyses/VoxelTuningAnalysis/slfmriGetVoxTuningParams.m) : wrapper to get voxel direction tuning parameters quantified by 1) calculating voxel tuning by taking the argmax response, 2) fitting a von Mises function to voxels responses or 3) calculating voxel directional vector.  
slfMRIwrapperPlotVoxTuningDB : plot database of voxel tuning (e.g., to stimulus motion direction) e.g., stimulus motion direction tuning curves for each voxel  
[slfmriGraphPermutedVoxTunings](fMRI/analyses/VoxelTuningAnalysis/slfmriGraphPermutedVoxTunings.m): visualize individual voxel tuning produced by permuted voxel responses.



##### voxel tuning population

[slfMRIwrapperDisVoxSelbalDir](fMRI/analyses/VoxelTuningAnalysis/slfMRIwrapperDisVoxSelbalDir.m) : plot the distribution of voxel selectivities for stimulus motion direction (a feature variable) sorted by another variable (the stimulus motion coherence, weak or strong), with the stimulus motion direction distributions equalized between the conditions of a third variable (switching : reported direction estimate near prior=1 or evidence=2).  

[slfMRIwrapperDisVoxSelbalDirVMfitfminS.m](fMRI/analyses/VoxelTuningAnalysis/slfMRIwrapperDisVoxSelbalDirVMfitfminS.m): get voxel tuning to motion direction by fitting a von mises function to voxel responses by motion direction using fminsearch to search for all the von mises parameters (amplitude, concentration and mean).

[slfmriGetVoxPopActivity](fMRI/analyses/VoxPopAnalysis/slfmriGetVoxPopActivity.m) : plot voxel population activity (response by voxel selectivities sorted by a variable, e.g., coherence)  


##### nested functions

[slfmrigetStimvolVars](MRI/analyses/slfmrigetStimvolVars.m) : get variable name and values associated with stimvols slfmriGetInstancedb : make a database of instances, corresponding variables and stimvols for one fMRI session  
[slmakeSwitchingVar](fMRI/analyses/classification/slmakeSwitchingVar.m) : classifies motion direction estimates into two classes (switched to prior and switched to direction) of a variable "mySwitch" and add that variable to the stimfiles associated with an fMRI session. The estimates is assigned to a class based on its vector distance to that class. e.g., Estimates nearest to statistical prior mean/motion direction are assigned value 1/2.  
[slfMRIwrapperDirDisBySwitch](fMRI/analyses/Controls/slfMRIwrapperDirDisBySwitch.m) : plot distribution of stimulus motion directions sorted by switching variable  

[slfMRIwrapperDirDisBySwitch_testSampleSize](fMRI/analyses/VoxelTuningAnalysis/slfMRIwrapperDirDisBySwitch_testSampleSize.m):  
[slfMRIwrapperDisVoxSelbalDir_testSampleSize](fMRI/analyses/VoxelTuningAnalysis/slfMRIwrapperDisVoxSelbalDir_testSampleSize.m):

##### fmri mapping  

SLfmriLocalizerMotion :
SLfmriProtocolanatRet_s03XX20150324  
SLfmriProtocolanatRetino_s030020150403  
SLfmriProtocolanatRetino_sFantom2015xxxx  
SLfmriProtocolanatRetinoMotionLoc_s025201501XX  
SLfmriProtocolanatRetinoMotionLoc_s02520150209  
SLfmriProtocolanatRetinoMotionLoc_s02520150224  
SLfmriProtocolanatRetinoMotionLoc_s02520150311  

##### Tasks 

SLtaskDotDirfMRI05: run the task during an fMRI scan for the "neural inference project"

