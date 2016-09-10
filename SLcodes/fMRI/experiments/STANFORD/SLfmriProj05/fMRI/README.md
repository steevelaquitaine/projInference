README

1 - Set experiment parameters

In matlab
```matlab
>> cd setParameters
>> edit SLinitRunExpUniPriorfMRI.m
```
Create a mat file containing the fMRI task parameters
e.g., here we draw 183 and 74 trials in which 5 stimulus motion directions are displayed with 6% and 12% interleaved motion coherences from a circular Gaussian (von Mises) distribution VM(theta:225,0.74848) with mean 225 degrees and concentration parameter k=0.74848 (about 80 degrees std)
```matlab
>> task = SLinitRunExpUniPriorfMRI('fMRIParamsPrior80Coh12and06',225,0.74848,[.06 .12],[183 74],15:70:355,[0.5 0 0])
```
