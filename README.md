# <img src="docs/mTRF-Toolbox_logo.png">

[![View mTRF-Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74260-mtrf-toolbox)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

mTRF-Toolbox is a MATLAB package for estimating neural encoding/decoding models and is suitable for neurophysiological data such as MEG, EEG, sEEG, ECoG and EMG. System identification is used to model how neuronal populations encode dynamic stimuli such as natural scenes and sounds. It can also be used to study other neural processes such as auditory/visual attention and multisensory integration, as well as various neural disorders where sensory processing is impaired. mTRF-Toolbox provides a basic machine learning framework for real-time BCI applications such as neural stimulus reconstruction and auditory attention decoding (AAD).

## Installation
Download and unzip mTRF-Toolbox to a local directory and in the MATLAB/GNU Octave command prompt enter:
```
addpath 'directory/mTRF-Toolbox-master'
savepath
```

## Documentation
For details on background and applications, please refer to the [mTRF-Toolbox paper](docs/Crosse_etal_FrontHumNeurosci_2016.pdf).

For examples of usage, please see [mTRF demos](mTRFdemos.m).

## mTRF Modelling Framework
Neural encoding of dynamic stimulus features can be quantified by computing forward encoding models, also known as temporal response functions or temporal receptive fields (TRFs). mTRF-Toolbox can process multivariate input features such as spatio- or spectro-temporal representations (STRFs), as well as categorical features such as phonetic or semantic embeddings. TRFs can be subjected to conventional time-frequency/source analysis techniques or used to predict the neural responses to new stimuli. Stimulus features can be reconstructed from neural responses by computing backward decoding models, whereby the direction of causality is reversed.
<div align="center">
  <img src="docs/mTRF-Toolbox.png" height="400">
</div>

## Contents
### Fitting Encoding/Decoding Models
* `mTRFcrossval()` - cross-validation for tuning model hyperparameters
* `mTRFtrain()` - encoding/decoding model fitting (TRF/STRF estimation)
* `mTRFpredict()` - model prediction and evaluation
* `mTRFtransform()` - transforms decoding models into neurophysiologically interpretable encoding models
 
### Decoding Attention and Multisensory Processing
* `mTRFattncrossval()` - cross-validation for building an attention decoder
* `mTRFmulticrossval()` - cross-validation for building an additive model of multisensory processing
* `mTRFmultitrain()` - additive multisensory model fitting (TRF/STRF estimation)

### Least Squares Estimation
* `olscovmat()` - ordinary least squares covariance matrix estimation
* `mlscovmat()` - multisensory least squares covariance matrix estimation

### Feature Extraction
* `mTRFenvelope()` - computes the acoustic envelope of an audio signal
* `mTRFresample()` - resamples and smooths temporal features
* `lagGen()` - generates time-lagged input features

## License
[BSD 3-Clause License](LICENSE)
