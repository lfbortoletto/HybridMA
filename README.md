HybridMA function corrects motion artifacts in time series data using a combination of spline interpolation and a filter based on the wavelet transform. The function identifies segments where motion artifacts exceed a specified threshold and applies a spline-based correction. Then, the spline corrected data is decomposed into wavelets, by wavelet transform, and coefficients are filtered based on wavelet threshold.

Inputs:

	timeSeriesData: Matrix of input data (time points x channels).
	acquisitionFrequency: Sampling frequency of the data acquisition system.
	splineThreshold: Motion detection identification threshold.
 	waveletThreshold: Parameter used to compute the statistics (iqr = 1.5 is 1.5 times the interquartile range and is usually used to detect outliers). Increasing it, it will delete fewer coefficients. If iqr<0 then this function is skipped. 
	K: Free parameter for motion correction (leave empty for default. Default: 2.5*acquisitionFrequency).
Outputs:

	hybridCorrectedData: Motion-corrected time series by hybrid approach, using spline and wavelet combination (time points x channels).

The algorithm follows the procedure found in,

Sergio L. Novi, Erin Roberts, Danielle Spagnuolo, Brianna M. Spilsbury, D'manda C. Price, Cara A. Imbalzano, Edwin Forero, Arjun G. Yodh, Glen M. Tellis, Cari M. Tellis, Rickson C. Mesquita, "Functional near-infrared spectroscopy for speech protocols: characterization of motion artifacts and guidelines for improving data analysis," Neurophoton. 7(1) 015001 (10 January 2020) https://doi.org/10.1117/1.NPh.7.1.015001

Please cite this work when using any of the codes.
