# ProteinClusterAnalysis

This is a workflow to for HS-AFM protein cluster dynamics analysis and in-cluster protein tilt and nearest-neighbor analysis. The codes are developed in BIO-AFM-LAB at Weill Cornell Medicine under the supervision of Professor Simon Scheuring.

Developer: Yining Jiang

Publication: XXXXX

User should email to the corresponding author of the paper for details about this work: Professor Simon Scheuring (sis2019@med.cornell.edu)

NOTE: Any usage of the codes should cite the publication mentioned above.

## System requirements:
1. Operating system for code development : macOS Big Sur Version 11.7.8
2. Software for code development: MATLAB (MathWorks) 2023b, Python 3.9.6
3. Additional add-ons: MIJI
4. Non-standard hardware: N/A

## Installation instructions: 
1. The codes require installation of MATLAB (MathWorks) 2023b. An installation guide can be found at: https://www.mathworks.com/help/install/.
2. MIJI is recommanded (not required) for visualizing data. An installation guide can be found at: https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab. If MIJI is not installed, user should comment out any code that uses MIJI for visualization (lines starting sith "MIJ.xxx").
3. The installation should take less than one hour on a "normal" desktop computer.

## General instructions:
These codes comprise of two parts:
### 1. HS-AFM protein clusters dynamics analysis
Note: These codes are developed to analyze the dynamic morphological changes of protein clusters imaged by HS-AFM
#### Main scripts:
1. analysis_radial_profile_main.m
2. analysis_radial_profile_core.m

#### Demo data
1. data_test_1.mat    This file contains test dataset (1) for a dynamic protein cluster
2. data_test_2.mat    This file contains test dataset (2) for a stable protein cluster 
3. data_test_result_1.mat    This file contains major outputs for test dataset 1
4. data_test_result_2.mat    This file contains major outputs for test dataset 2

#### Instruction
User should run the main scripts. Operation details are provided within the scripts.

### 2. HS-AFM protein tilt and nearest-neighbor analysis
Note: These codes are developed to analyze the dynamic morphological changes of protein clusters imaged by HS-AFM
#### Main scripts:
1. patch_tilt_nn_analysis.m

#### Demo data
1. data_test.mat    This file contains test dataset for a cluster of 7 pentamer proteins
2. data_test_result.mat    This file contains major outputs for the test dataset

#### Instruction
User should run the main scripts. Operation details are provided within the scripts.


## Demo (Test data)
Test data for each part is provided as a 'xxx.mat' file. Essential output files could be found in 'xxx_result.mat' files.
These tests are expected to run for less than 5 minutes for demo on a "normal" desktop computer following the instructions provided in the main scripts.

## Instruction for use
User should process raw HS-AFM movies, including proper flattening and denoising. For HS-AFM protein clusters dynamics analysis, user should apply mask the pixels corresponding to the protein cluster (which should be the major feature in the movie). For HS-AFM protein tilt and nearest-neighbor analysis (of stable protein cluster), user should average the HS-AFM movie stack to obtain a single frame of the protein cluster. Besides, user should provide the coordinates of the single-particles and the molecular symmetry values. Details can be found alongside the scripts accordingly. User should adjust the parameters in the main scripts accourding to the nature of their proteins of interest.
