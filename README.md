[![View mTRF-Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74260-mtrf-toolbox)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# mTRF-Toolbox
mTRF-Toolbox is a MATLAB package for rapid estimation of forward encoding
models (stimulus to neural response) or backward decoding models (neural
response to stimulus), suitable for modelling neurophysiological data such 
as MEG, EEG, iEEG, sEEG, ECoG and EMG data. 

Forward encoding models, also known as response functions or 
receptive fields, can be used to investigate information processing in 
neuronal populations with respect to temporal features (TRFs), or 
spectro- or spatio-temporal features (STRFs). STRFs can be subjected to 
conventional time-frequency and source analysis techniques used to analyse
event related potentials. In addition, TRFs can be used to predict
the dynamics of neural responses to unseen stimuli as a way to objectively 
measure stimulus encoding. Stimulus reconstruction can be performed using 
backward decoding models that project the multi-channel neural responses 
back to the dynamics of the stimulus. This is useful for decoding stimulus 
features from neural responses and can be used to build brain-computer 
interfaces and other real-time neural interface applications.

mTRF-Toolbox facilitates the use of natural continuous stimuli, allowing 
researchers to investigate how neural systems process dynamic environmental 
signals such as speech, music and motion, and to decode dynamic cognitive 
processes such as auditory attention and multisensory integration.

### Documentation
Documentation on mTRF-Toolbox usage and underlying theory can be found [here](http://mickcrosse.com/assets/pubs/Crosse_etal_FrontHumNeurosci_2016.pdf).

## mTRF Modelling Framework
<img src="doc/mTRF-Toolbox.png" width="600" height="400">

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

### Rapid Least Squares Estimation
* `olscovmat()` - ordinary least squares covariance matrix estimation
* `mlscovmat()` - multisensory least squares covariance matrix estimation

### Feature Extraction
* `mTRFenvelope()` - computes the acoustic envelope of an audio signal
* `mTRFresample()` - resamples and smooths temporal features
* `lagGen()` - generates time-lagged input features

## License
[BSD 3-Clause License](LICENSE)
