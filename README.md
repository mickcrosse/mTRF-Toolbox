# <img src="doc/mTRF-Toolbox_logo.png">

[![View mTRF-Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74260-mtrf-toolbox)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

mTRF-Toolbox is a MATLAB package for estimating neural encoding/decoding models, suitable for neurophysiological data such as MEG, EEG, sEEG, ECoG and EMG.

Neural encoding of dynamic stimulus features can be quantified by computing forward encoding models, also known as temporal response functions or receptive fields (TRFs). mTRF-Toolbox can process multivariate input features such as spatio- or spectro-temporal representations (STRFs), as well as categorical features such as phonemes or semantics. TRFs can be subjected to conventional time-frequency/source analysis techniques or used to predict the neural responses to new stimuli. Stimulus features can be reconstructed from the neural responses by computing backward decoding models, whereby the direction of causality is reversed.

With mTRF-Toolbox, it is possible to study how neuronal populations encode natural scenes and sounds such as motion, speech and music. It can also be used to study neural processes such as auditory/visual attention and multisensory integration, as well as neural disorders that impair sensory processing. The mTRF modelling framework provides the basic tools for building brain-computer interfaces and other real-time neural interface applications.

### Documentation
Documentation on toolbox usage and theory is available in the [mTRF-Toolbox paper](http://mickcrosse.com/assets/pubs/Crosse_etal_FrontHumNeurosci_2016.pdf).

## mTRF Modelling Framework
<img src="doc/mTRF-Toolbox.png" height="400">

## Contents
### Fitting Encoding/Decoding Models
* `mTRFcrossval()` - cross-validation for tuning model hyperparameters
* `mTRFtrain()` - encoding/decoding model fitting (TRF/STRF estimation)
* `mTRFpredict()` - model prediction and evaluation
* `mTRFtransform()` - transforms decoding models into neurophysiologically interpretable encoding models
 
### Decoding Attention and Multisensory Processing
* `mTRFattncrossval()` - cross-validation for bulding an attention decoder
* `mTRFmulticrossval()` - cross-validation for bulding an additive model of multisensory processing
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
