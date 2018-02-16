mTRF Toolbox is a MATLAB TOOLBOX that permits the fast computation of the
linear stimulus-response mapping of any sensory system in the forward or 
backward direction. It is suitable for analysing EEG, MEG, ECoG and EMG 
data.

The forward model, or temporal response function (TRF), can be interpreted
using conventional analysis techniques such as time-frequency and source 
analysis. The TRF can also be used to predict future responses of the 
system given a new stimulus signal. Similarly, the backward model can be 
used to reconstruct spectrotemporal stimulus information given new response 
data.

mTRF Toolbox facilitates the use of continuous stimuli in
electrophysiological studies as opposed to time-locked averaging
techniques which require discrete stimuli. This enables examination of
how neural systems process more natural and ecologically valid stimuli
such as speech, music, motion and contrast.

Table of Contents
=================

* mTRFtrain.m Usage
* mTRFpredict.m Usage
* mTRFtransform.m Usage
* mTRFcrossval.m Usage
* mTRFmulticrossval.m Usage
* lagGen.m Usage
* Sample Data Sets
* Tips on practical Use
* Examples
* Additional Information

mTRFtrain Usage
===============

mTRFtrain mTRF Toolbox training function.
  MODEL = MTRFTRAIN(STIM,RESP,FS,MAP,TMIN,TMAX,LAMBDA) performs ridge
  regression on the stimulus property STIM and the neural response data
  RESP to solve for their linear mapping function MODEL. Pass in MAP==1
  to map in the forward direction or MAP==-1 to map backwards. The
  sampling frequency FS should be defined in Hertz and the time lags
  should be set in milliseconds between TMIN and TMAX. Regularisation is
  controlled by the ridge parameter LAMBDA.

  [...,T,C] = MTRFTRAIN(...) also returns the vector of time lags T for
  plotting MODEL and the regression constant C for absorbing any bias
  when testing MODEL.

  Inputs:
  stim   - stimulus property (time by features)
  resp   - neural response data (time by channels)
  fs     - sampling frequency (Hz)
  map    - mapping direction (forward==1, backward==-1)
  tmin   - minimum time lag (ms)
  tmax   - maximum time lag (ms)
  lambda - ridge parameter

  Outputs:
  model  - linear mapping function (MAP==1: feats by lags by chans,
           MAP==-1: chans by lags by feats)
  t      - vector of time lags used (ms)
  c      - regression constant

mTRFpredict Usage
=================

mTRFpredict mTRF Toolbox prediction function.
  PRED = MTRFPREDICT(STIM,RESP,MODEL,FS,MAP,TMIN,TMAX,C) performs a
  convolution of the stimulus property STIM or the neural response data
  RESP with their linear mapping function MODEL to solve for the
  prediction PRED. Pass in MAP==1 to predict RESP or MAP==-1 to predict
  STIM. The sampling frequency FS should be defined in Hertz and the time
  lags should be set in milliseconds between TMIN and TMAX. The
  regression constant C absorbs any bias in MODEL.

  [...,R,P,MSE] = MTRFPREDICT(...) also returns the correlation
  coefficients R between the original and predicted values, the
  corresponding p-values P and the mean squared errors MSE.

  Inputs:
  stim   - stimulus property (time by features)
  resp   - neural response data (time by channels)
  model  - linear mapping function (MAP==1: feats by lags by chans,
           MAP==-1: chans by lags by feats)
  fs     - sampling frequency (Hz)
  map    - mapping direction (forward==1, backward==-1)
  tmin   - minimum time lag (ms)
  tmax   - maximum time lag (ms)
  c      - regression constant

  Outputs:
  pred   - prediction (MAP==1: time by chans, MAP==-1: time by feats)
  r      - correlation coefficients
  p      - p-values of the correlations
  mse    - mean squared errors

mTRFtransform Usage
===================

mTRFtransform mTRF Toolbox mapping transformation function.
   MODEL_T = MTRFTRANSFORM(B_MODEL,RESP,STIM) tansforms the coefficients 
   of the model weights MODEL into transformed model coefficients MODEL_T.
 
   Inputs:
   stim   - stimulus property (time by features)
   resp   - neural response data (time by channels)
   model  - linear mapping function (MAP==1: feats by lags by chans,
            MAP==-1: chans by lags by feats)
   fs     - sampling frequency (Hz)
   map    - transformation direction (forward -> backward==1, backward -> 
            forward==-1)
   tmin   - minimum time lag (ms)
   tmax   - maximum time lag (ms)
   c      - regression constant

   Outputs:
   model_t - transformed model weights (lags by chans)
   t       - vector of time lags used (ms)
   c_t     - transformed model constant

mTRFcrossval Usage
==================

mTRFcrossval mTRF Toolbox cross-validation function.
  [R,P,MSE] = MTRFCROSSVAL(STIM,RESP,FS,MAP,TMIN,TMAX,LAMBDA) performs
  leave-one-out cross-validation on the set of stimuli STIM and the
  neural responses RESP for the range of ridge parameter values LAMBDA.
  As a measure of performance, it returns the correlation coefficients R
  between the predicted and original signals, the corresponding p-values
  P and the mean squared errors MSE. Pass in MAP==1 to map in the forward
  direction or MAP==-1 to map backwards. The sampling frequency FS should
  be defined in Hertz and the time lags should be set in milliseconds
  between TMIN and TMAX.

  [...,PRED,MODEL] = MTRFCROSSVAL(...) also returns the predictions PRED
  and the linear mapping functions MODEL.

  Inputs:
  stim   - set of stimuli [cell{1,trials}(time by features)]
  resp   - set of neural responses [cell{1,trials}(time by channels)]
  fs     - sampling frequency (Hz)
  map    - mapping direction (forward==1, backward==-1)
  tmin   - minimum time lag (ms)
  tmax   - maximum time lag (ms)
  lambda - ridge parameter values

  Outputs:
  r      - correlation coefficients
  p      - p-values of the correlations
  mse    - mean squared errors
  pred   - prediction [MAP==1: cell{1,trials}(lambdas by time by chans),
           MAP==-1: cell{1,trials}(lambdas by time by feats)]
  model  - linear mapping function (MAP==1: trials by lambdas by feats by
           lags by chans, MAP==-1: trials by lambdas by chans by lags by
           feats)

mTRFmulticrossval Usage
=======================

mTRFmulticrossval mTRF Toolbox multisensory cross-validation function.
  [R,P,MSE] = MTRFMULTICROSSVAL(STIM,RESP,RESP1,RESP2,FS,MAP,TMIN,TMAX,
  LAMBDA1,LAMBDA2) performs leave-one-out cross-validation of an
  additive model for a multisensory dataset as follows:
  1. Separate unisensory models are calculated using the set of stimuli
     STIM and unisensory neural responses RESP1 and RESP2 for the range
     of ridge parameter values LAMBDA1 and LAMBDA2 respectively.
  2. The algebraic sums of the unisensory models for every combination of
     LAMBDA1 and LAMBDA2 are calculated, i.e., the additive models.
  3. The additive models are validated by testing them on the set of
     multisensory neural responses RESP.
  As a measure of performance, it returns the correlation coefficients R
  between the predicted and original signals, the corresponding p-values
  P and the mean squared errors MSE. The time lags T should be set in
  milliseconds between TMIN and TMAX and the sampling frequency FS should
  be defined in Hertz. Pass in MAP==1 to map in the forward direction or
  MAP==-1 to map backwards. The neural responses in all three sensory
  conditions must have been recorded for the same set of stimuli STIM.

  [...,PRED,MODEL] = MTRFMULTICROSSVAL(...) also returns the predictions
  PRED and the linear mapping functions MODEL.

  Inputs:
  stim    - set of stimuli [cell{1,trials}(time by features)]
  resp    - set of multisensory neural responses [cell{1,trials}(time by channels)]
  resp1   - set of unisensory 1 neural responses [cell{1,trials}(time by channels)]
  resp2   - set of unisensory 2 neural responses [cell{1,trials}(time by channels)]
  fs      - sampling frequency (Hz)
  map     - mapping direction (forward==1, backward==-1)
  tmin    - minimum time lag (ms)
  tmax    - maximum time lag (ms)
  lambda1 - unisensory 1 ridge parameter values
  lambda2 - unisensory 2 ridge parameter values

  Outputs:
  r       - correlation coefficients
  p       - p-values of the correlations
  mse     - mean squared errors
  pred    - prediction [MAP==1: cell{1,trials}(lambdas1 by lambdas2 by
            time by chans), MAP==-1: cell{1,trials}(lambdas1 by lambdas2
            by time by feats)]
  model   - linear mapping function (MAP==1: trials by lambdas1 by
            lambdas2 by feats by lags by chans, MAP==-1: trials by
            lambdas1 by lambdas2 by chans by lags by feats)

lagGen Usage
============

lagGen Lag generator.
  [XLAG] = LAGGEN(X,LAGS) returns the matrix XLAG containing the lagged
  time series of X for a range of time lags given by the vector LAGS. If
  X is multivariate, LAGGEN will concatenate the features for each lag
  along the columns of XLAG.

  Inputs:
  x    - vector or matrix of time series data (time by features)
  lags - vector of integer time lags (samples)

  Outputs:
  xLag - matrix of lagged time series data (time by lags*feats)

Sample Data Sets
================

contrast_data.mat
This MATLAB file contains 3 variables. The first is a matrix consisting
of 120 seconds of 128-channel EEG data. The second is a vector consisting
of a normalised sequence of numbers that indicate the contrast of a
checkerboard that was presented during the EEG at a rate of 60 Hz. The
third is a scaler which represents the sample rate of the contrast signal
and EEG data (128 Hz). See Lalor et al. (2006) for further details.

coherentMotion_data.mat
This MATLAB file contains 3 variables. The first is a matrix consisting
of 200 seconds of 128-channel EEG data. The second is a vector consisting
of a normalised sequence of numbers that indicate the motion coherence of
a dot field that was presented during the EEG at a rate of 60 Hz. The
third is a scaler which represents the sample rate of the motion signal
and EEG data (128 Hz). See Gonçalves et al. (2014) for further details.

speech_data.mat
This MATLAB file contains 4 variables. The first is a matrix consisting
of 120 seconds of 128-channel EEG data. The second is a matrix consisting
of a speech spectrogram. This was calculated by band-pass filtering the 
speech signal into 128 logarithmically-spaced frequency bands between 100 
and 4000 Hz and taking the Hilbert transform at each frequency band. The 
spectrogram was then downsampled to 16 frequency bands by averaging 
across every 8 neighbouring frequency bands. The third variable is the 
broadband envelope, obtained by taking the mean across the 16 narrowband 
envelopes. The fourth variable is a scaler which represents the sample 
rate of the envelope, spectrogram and EEG data (128 Hz). See Lalor & 
Foxe (2010) for further details.

Tips on Practical Use
=====================

* Ensure that the stimulus and response data have the same sample rate
  and number of samples.
* Downsample the data when conducting large-scale multivariate analyses
  to reduce running time, e.g., 128 Hz or 64 Hz.
* Normalise all data, e.g., between [-1,1] or [0,1] or z-score. This will 
  stabalise regularisation across trials and enable a smaller parameter 
  search.
* Enter the start and finish time lags in milliseconds. Enter positive
  lags for post-stimulus mapping and negative lags for pre-stimulus
  mapping. This is the same for both forward and backward mapping - the 
  code will automatically reverse the lags for backward mapping.
* When using MTRFPREDICT, always enter the model in its original
  3-dimensional form, i.e., do not remove any singleton dimensions.
* When using MTRFCROSSVAL, the trials do not have to be the same length,
  but using trials of the same length will optimise performance.
* When using mTRFmulticrossval, the trials in each of the three sensory
  conditions should correspond to the stimuli in STIM.

Examples
========

Contrast: Forward model (TRF/VESPA)
>> load('contrast_data.mat');
>> [w,t] = mTRFtrain(contrastLevel,EEG,128,0,-150,450,1);
>> figure; plot(t,squeeze(w(1,:,23))); xlim([-100,400]);
>> xlabel('Time lag (ms)'); ylabel('Amplitude (a.u.)')

Motion: Forward model (TRF/VESPA)
>> load('coherentmotion_data.mat');
>> [w,t] = mTRFtrain(coherentMotionLevel,EEG,128,0,-150,450,1);
>> figure; plot(t,squeeze(w(1,:,21))); xlim([-100,400]);
>> xlabel('Time lag (ms)'); ylabel('Amplitude (a.u.)')

Speech: Forward model (TRF/AESPA)
>> load('speech_data.mat');
>> [w,t] = mTRFtrain(envelope,EEG,128,0,-150,450,0.1);
>> figure; plot(t,squeeze(w(1,:,85))); xlim([-100,400]);
>> xlabel('Time lag (ms)'); ylabel('Amplitude (a.u.)')

Speech: Spectrotemporal forward model (STRF)
>> load('speech_data.mat');
>> [w,t] = mTRFtrain(spectrogram,EEG,128,0,-150,450,100);
>> figure; imagesc(t,1:16,squeeze(w(:,:,85))); xlim([-100,400]);
>> xlabel('Time lag (ms)'); ylabel('Frequency band')

Speech: Backward model (stimulus reconstruction)
>> load('speech_data.mat');
>> envelope = resample(envelope,64,128); EEG = resample(EEG,64,128);
>> stimTrain = envelope(1:64*60,1); respTrain = EEG(1:64*60,:);
>> [g,t,con] = mTRFtrain(stimTrain,respTrain,64,1,0,500,1e5);
>> stimTest = envelope(64*60+1:end,1); respTest = EEG(64*60+1:end,:);
>> [recon,r,p,MSE] = mTRFpredict(stimTest,respTest,g,64,1,0,500,con);

Additional Information
======================

mTRF Toolbox is available for download at:
http://sourceforge.net/projects/aespa

mTRF Toolbox support documentation is available at:
http://dx.doi.org/10.3389/fnhum.2016.00604

For any questions and comments, please email Dr. Edmund Lalor at:
edmundlalor@gmail.com

Acknowledgments:
This work was supported in part by Science Foundation Ireland (SFI), the
Irish Higher Education Authority (HEA) and the Irish Research Council
(IRC).

References:
* Lalor EC, Pearlmutter BA, Reilly RB, McDarby G, Foxe JJ (2006) The
  VESPA: a method for the rapid estimation of a visual evoked potential.
  NeuroImage 32:1549-1561.
* Gonçalves NR, Whelan R, Foxe JJ, Lalor EC (2014) Towards obtaining
  spatiotemporally precise responses to continuous sensory stimuli in
  humans: a general linear modeling approach to EEG. NeuroImage 97(2014):196-205.
* Lalor, EC, & Foxe, JJ (2010) Neural responses to uninterrupted natural
  speech can be extracted with precise temporal resolution. Eur J Neurosci 
  31(1):189-193.
* Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2015) The multivariate 
  temporal response function (mTRF) toolbox: a MATLAB toolbox for relating 
  neural signals to continuous stimuli. Front Hum Neurosci 10:604.
