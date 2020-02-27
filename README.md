# <img src="docs/mTRF-Toolbox_logo.png">

[![View mTRF-Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74260-mtrf-toolbox)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

mTRF-Toolbox is a MATLAB package for estimating linear neural encoding/decoding models, suitable for neurophysiological data such as MEG, EEG, sEEG, ECoG and EMG. System identification is used to model how neuronal populations encode dynamic stimuli such as natural scenes and sounds. It can also be used to study other neural processes such as auditory/visual attention and multisensory integration, as well as various neural disorders where sensory processing is impaired. mTRF-Toolbox provides a basic machine learning framework for real-time BCI applications such as neural stimulus reconstruction and auditory attention decoding (AAD).

- [Installation](#installation)
- [Documentation](#documentation)
- [mTRF Modelling Framework](#mtrf-modelling-framework)
- [Contents](#contents)
  - [Fitting Encoding and Decoding Models](#fitting-encoding-and-decoding-models)
  - [Decoding Attention and Multisensory Processing](#decoding-attention-and-multisensory-processing)
  - [Least Squares Estimation](#least-squares-estimation)
  - [Feature Extraction](#feature-extraction)
- [License](#license)

## Installation
Download and unzip mTRF-Toolbox to a local directory, then in the MATLAB/GNU Octave command prompt enter:
```
addpath 'directory/mTRF-Toolbox-master'
savepath
```

## Documentation
For documentation, please refer to the [mTRF-Toolbox paper](docs/Crosse_etal_FrontHumNeurosci_2016.pdf).

For code demonstrating usage, please see [mTRFdemos](mTRFdemos.m).

## mTRF Modelling Framework
mTRF-Toolbox provides a complementary forward (encoding) and backward (decoding) modelling framework. Known as a temporal response function or temporal receptive field (TRF), a forward model describes the process by which sensory information is encoded in the neural activity. Multivariate stimulus features such as spatio- or spectro-temporal representations, as well as categorical features such as phonetic or semantic embeddings can be used as inputs to the model. TRFs can be subjected to conventional time-frequency/source analysis techniques or used to predict the neural responses to an independent set of stimuli. Known as a decoder or reconstruction filter, a backward model reverses the direction of causality between stimulus and response. mTRF-Toolbox provides an efficient cross-validation framework for optimizing hyperparameter configuration. Backward models can be used to reconstruct stimulus features from an independent set of neural responses, or decode higher-order cognitive processes such as top-down attention. 

<div align="center">
  <img src="docs/mTRF-Toolbox.png" height="400">
</div>

## Contents
### Fitting Encoding and Decoding Models
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
