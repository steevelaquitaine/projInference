## README

 * author: Steeve Laquitaine
 * date: 160424
 * purpose:  run motion direction estimation task in fMRI where the statistical distribution of motion directions and the coherence of the motion change.
 * Description: Subjects are first trained with one prior distribution then scanned with the same prior. In subsequent session, subjects are trained with a new prior and scanned with this new prior. 
 * Programming language: Matlab

To run the experiment step-by-step (in matlab)

```matlab
  >> edit runAllTaskDirfMRI05.m
```

#### Example Experiments

Probe motion direction statistics induced-changes in voxel direction selectivity distributions

1 direction statistics (V(directions;u=225,k=0.74848)) run in three scans and 2 weak motions (low coherences) interleaved within scans

sub01 - id: s025 : 3 scans, coh:6%,12%, : prior mean=225 (150814; 150923; 150925); 1 scan, coh:6%,12%, : prior mean=135 (160614)

Results: We can use machine learning to decode (classify) coherence, direction and when subject switch to the prior or to the direction.



