This function corrects motion artifacts in time series data using spline interpolation. The function identifies segments where motion artifacts exceed a specified threshold and applies a spline-based correction.

Inputs:
	timeSeriesData: Matrix of input data (time points x channels)
	acquisitionFrequency: Sampling frequency of the data acquisition system.
	threshold: Motion detection threshold.
	K: Free parameter for motion correction (leave empty for default. Default: 2.5*acquisitionFrequency).
Outputs:
	correctedData: Motion-corrected time series.

The algorithm follows the procedure found in,

Sergio L. Novi, Erin Roberts, Danielle Spagnuolo, Brianna M. Spilsbury, D'manda C. Price, Cara A. Imbalzano, Edwin Forero, Arjun G. Yodh, Glen M. Tellis, Cari M. Tellis, Rickson C. Mesquita, "Functional near-infrared spectroscopy for speech protocols: characterization of motion artifacts and guidelines for improving data analysis," Neurophoton. 7(1) 015001 (10 January 2020) https://doi.org/10.1117/1.NPh.7.1.015001

Please cite this work when using any of the codes.
