# PVF Vasomotion Analysis

This repository contains MATLAB analysis scripts used to quantify perivascular fibroblast (PVF) calcium signaling and pial artery vasomotor dynamics for the study **"Perivascular fibroblasts locally regulate the vasomotor dynamics of pial arterioles,**" performed by the Bonney Lab at the University of Colorado Anschutz Medical Campus.

These scripts were developed for analysis of two-photon _in vivo_ imaging data examining relationships between cellular calcium activity and vascular dynamics.

## About the Study

Cerebral arteries and arterioles undergo spontaneous rhythmic changes in vessel diameter known as vasomotion. In addition to contributing to the regulation of cerebral blood flow, these vascular dynamics are associated with cerebrospinal fluid (CSF) movement and perivascular waste clearance.

This study investigates the relationship between perivascular fibroblast calcium signaling and pial artery vasomotor dynamics, including the temporal and phase relationships between cellular calcium activity and vessel diameter oscillations, synchronization of vasomotor oscillations along the length of arterioles, and vasomotor frequency metrics.

## Publication

The code in this repository is associated with the following publication:

**[citation here].**

DOI: **[DOI here]**

## Repository Overview

The repository is organized into four analysis categories:

### Calcium_Vasomotion

Analyses examining cellular calcium dynamics in relation to vasomotor oscillatory dynamics.

- `Calcium_Vaso_Amplitude.m`  
  Calculates maximum, minimum, and amplitude values for calcium and vessel diameter signals.

- `Calcium_Vaso_CrossCorrelation.m`  
  Performs cross-correlation analysis between calcium and vessel diameter signals.

- `Calcium_Vaso_CycleTriggeredAverage.m`  
  Calculates cycle-triggered average of calcium dynamics aligned to vasomotor oscillations.

- `Calcium_Vaso_FrequencyAnalysis.m`  
  Calculates calcium and vessel diameter frequency metrics using Welch power spectral density and Hilbert-based instantaneous frequency analyses.

- `Calcium_Vaso_HilbertPhase.m`  
  Calculates the phase relationship between calcium and vessel diameter oscillations using Hilbert transform-based phase analysis.

### Vessel_Synchrony

Analyses examining synchronization of vasomotor oscillations between ROIs along the length of vessels.

- `Vessel_Coherence.m`  
  Calculates pairwise magnitude-squared coherence between vessel diameter ROIs.

- `VesselROI_FrequencyAnalysis.m`  
  Calculates frequency metrics of individual vessel diameter ROIs using Welch power spectral density and Hilbert-based instantaneous frequency analyses.

- `VesselSynchrony_Hilbert.m`  
  Calculates pairwise phase synchrony between vessel diameter ROIs using Hilbert phase analysis and phase-locking value (PLV).

### Vessel_Diameter

- `VesselDiameter_FrequencyAnalysis.m`  
  Calculates vessel diameter frequency characteristics using Welch power spectral density and Hilbert-based instantaneous frequency analyses.

### Circular_Statistics

- `CircularStatTests_WW_MWW.m`  
  Performs Watson-Williams and Mardia-Watson-Wheeler tests for comparisons of circular data across groups, including pairwise comparisons with Holm correction for multiple comparisons.

## Requirements

These analyses were developed in MATLAB.

Required MathWorks toolboxes:

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

MATLAB version used for script development and analyses: **[R2025b]**

## Input Data

Most analyses use `.csv` files in a standardized format where row 1 contains column headers, rows 2–11 contain acquisition metadata from two-photon imaging, and numeric time-series data begin on row 12.

For calcium and vessel diameter analyses, the standard input format is:

- Column E: time
- Column F: calcium signal (ΔF/F0)
- Column G: vessel diameter (% baseline)

For vessel synchrony analyses:

- Column E: time
- Column F onward: vessel diameter ROIs (% baseline)

For circular statistics analysis:

- (`.xlsx`) files in which each column represents a group, with individual circular mean values from Hilbert analyses expressed in degrees.

Input files should be given descriptive names to help identify and organize analysis outputs when batch processing, for example:

`Animal1_Area1_ROI1.csv`

Example input files demonstrating the standardized format are provided in the [`Example_Input`] folder.

- [`calcium_vasomotion_example_input.csv`] – for calcium and vessel diameter analyses.

- [`vessel_synchrony_example_input.csv`] – for vessel synchrony analyses.

- [`circular_statistics_example_input.xlsx`] – for circular statistical analyses.

## Setup Instructions

### 1. Download the repository

Download the repository using **Code → Download ZIP** and extract the files, or clone the repository using Git.

Then, open the analysis scripts in MATLAB.

### 2. Specify the input directory

Before running any script, replace the placeholder input directory (shown below) near the top of the script with the path to the folder containing the input data:

```matlab
inputFolder = 'input\directory\here'; % Replace with path to user input data
```

The scripts will recursively search for all .csv files in the specified input directory and its subfolders for input files.

The circular statistics script searches the specified input directory for Excel (`.xlsx`) files.

### 3. Check analysis settings

Analysis parameters are defined near the top of each script. Confirm that the settings are appropriate for the dataset before running the analysis.

Depending on the script, adjustable parameters could include:

- Sampling frequency (`Fs`)
- Frequency range
- Filter settings
- Detrending
- Edge trimming
- Peak detection settings
- Cross-correlation lag window
- Welch PSD settings
- Hilbert amplitude threshold

Sampling frequency should be set according to the acquisition rate of the input recording. Sampling frequencies used across these analyses include approximately 1.7 Hz, 2 Hz, and 7 Hz recordings. The specific values used for each analysis are commented within the scripts on the line corresponding to sampling frequency settings, which will appear as seen below:

```matlab
Fs = 1.951; % Sampling frequency: 1.951 for 2Hz, 1.727 for 1.7Hz, 7.440476 for 7Hz
```

### 4. Run the analysis

Run the desired `.m` file in MATLAB.

Scripts process the input files and generate analysis-specific output files and folders. Output locations are defined within each script and will generally be found in the chosen input directory.

## Analysis Outputs

Output files are saved as `.csv`, `.xlsx`, and/or figure files depending on the analysis.

## Citation

If you use these scripts or data in your research, please cite the original associated publication:

**[citation here]**

## Authors

**Stephanie Karas**  
Bonney Lab  
University of Colorado Anschutz Medical Campus

Additional contributors: **Bonney Lab: (https://github.com/BonneyLabatCU)**

## Contact

For inquiries or issues, please contact:

**Stephanie Karas**  
**stephanie.karas@cuanschutz.edu**
