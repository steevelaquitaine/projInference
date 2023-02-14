# A Switching Bayesian observer

Steeve Laquitaine & Justin L. Gardner

<p align="center">
	<img src="assets/switching_observer.png" width="800">
</p>

Analyse and model human inference behavior in motion direction and spatial location estimation tasks in which optimality requires exploiting knowledge of the statistical distribution of the stimuli.

See also [website](http://steevelaquitaine.github.io/projInference/).

## Task

To run the task git clone the project and follow the instructions given in .../task/README.txt

1. Set your screen parameters (will run on screenNumber = 1, you can change the screen number in taskDotDir.m)

2. Open runTask
	
Run an  10º prior

```matlab
taskDotDir('steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_073_033perCoh_130217')
```

or a 80º prior

```matlab
taskDotDir('steeve_exp12_metho_Pstd080_mean225_coh006012024_dir36_t106_075_034perCoh_130217')%dropped frames: 0.00 %done
```

note : taskDotDir.m loads the stimulus parameters (each trial's motion direction and coherence and a few stats) stored in a file like (e.g., 10º prior) :

```
steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_074_033perCoh_controlRndInitPosSymPrior_131224.mat
```

There is one file per prior.

3. When the task is running to enter a response move your mouse to the desired direction, then press 1 on keyboard (make sure mglEditScreenParameters takes "1" as input)


## Data extraction

1. Download the dataset from Mendeley into your project path e.g., `proj_path = "Desktop/project"`

```python
Desktop/project/
- data/
  - data01_direction4priors/ # four-prior motion experiment
      - data/
        - sub01/  # subject 1
          - steeve_exp12_data_sub01_sess06_run26_Pstd040_mean225_coh006_012_024_dir36_randInitPosSymPrior80_140106.mat # task & behavioral data
          - steeve_exp12_data_sub01_sess06_run26_Pstd040_mean225_coh006_012_024_dir36_randInitPosSymPrior80_140106.edf # eye tracking data
          ...
        - sub12/
  - data02_direction1prior/ # one-prior motion experiment
      - data/
        - sub01/
          - steeve_data_sub01_sess01_run01_Pstd080_mean225_coh006_012_024_dir36_170704.mat
          - steeve_data_sub01_sess01_run02_Pstd080_mean225_coh006_012_024_dir36_170704.edf
          ...
        - sub06/
  - data03_orientation/ # four-prior location experiment
      - data/
        - sub01/
          - steeve_exp04_data_sub01_sess01_run01_Pstd080_mean225_con010_0156_1_loc36perCon_151208.mat 
          - steeve_exp04_data_sub01_sess01_run01_Pstd080_mean225_con010_0156_1_loc36perCon_151208.edf
          ...
        - sub09/  
```

2. Extract the task variables from the data files:

The dataset is saved in `stim_file.mat`. You can access its variables with `getTaskParameters('stim_file.mat')`:

Clone the `mgl` software in your project path. In Matlab:

```python
# move to project path
cd Desktop/project

# clone package
git clone https://github.com/justingardner/mgl.git mgl

# add package to path
addpath(genpath('Desktop/project/mgl'))

# test your installation
help mgl
# MGL library functions
#
# Main screen functions
#   mglOpen                   : Opens the screen
#   mglFlush                  : Flips front and back buffer
#   mglClose                  : Closes the screen
```

'/gpfs/bbp.cscs.ch/data/scratch/proj68/laquitai/other_datasets/switching_data/data01_direction4priors/data/sub01'

## References 

Please cite:

```
@article{laquitaine2018switching,
  title={A switching observer for human perceptual estimation},
  author={Laquitaine, Steeve and Gardner, Justin L},
  journal={Neuron},
  volume={97},
  number={2},
  pages={462--474},
  year={2018},
  publisher={Elsevier}
}
```
