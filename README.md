# kumar-sensorimotor-integration

## This is the code used to analyze the data presented in S. Kumar, A. K. Sharma, A. Tran, and A. M. Leifer, titled "Inhibitory feedback from the motor circuit gates mechanosensory processing in *Caenorhabditis elegans*". 

The code is roughly divided into 2 sections: Real-time recording and optogenetic stimulation in LabVIEW, and analyzing data to generate figures in MATLAB

### Section 1: LabVIEW Real-time

This study uses a targeted optogenetic delivery system published previously. More details can be found at Liu et al., *PLoS Biology* 2022 (https://doi.org/10.1371/journal.pbio.3001524).

Requirements:
 - LabVIEW 2019 64-bit Windows
 - LabVIEW package manager
 - LabVIEW OpenG add-on
 - HDF5 v1.8.18+ (latest v1.8 but not v1.10)
 - h5labview

To run the real-time LabVIEW code, open the LabVIEW project at `\LabviewVIs\ProjectAPI.lvproj` and then select the appropriate experimental protocol vi. For example, open loop stimulation is done using `RunFullWormRails.vi`, closed-loop stimulation on turns was done using `RunRailsTriggeredByTurning.vi`, closed-loop stimulation on reversals was done using `RunStimulateReversingWorms.vi`. 

### Section 2: Generating figures

The MATLAB script `analysis_code_figure_1_2_4_S1_S3.m` is used to generate figures 1, 2, 4, S1, and S3. The script `analysis_code_figure_3_S2.m` generates figures 3 and S2 of the manuscript.