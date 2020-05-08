# <img src="docs/mTRF-Toolbox_logo.png">

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![View mTRF-Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74260-mtrf-toolbox)
[![Download mTRF-Toolbox](https://img.shields.io/sourceforge/dt/aespa.svg)](https://sourceforge.net/projects/aespa/files/latest/download)

mTRF-Toolbox is a MATLAB package for quantitative modelling of sensory processing, suitable for neurophysiological data such as MEG, EEG, sEEG, ECoG and EMG. It can be used to model the functional relationship between neuronal populations and dynamic sensory inputs such as natural scenes and sounds, or build neural decoders for reconstructing stimulus features and developing real-time applications such as brain-computer interfaces (BCIs).

- [Installation](#installation)
- [Documentation](#documentation)
- [mTRF Modelling Framework](#mtrf-modelling-framework)
- [Contents](#contents)
  - [Fitting encoding/decoding models](#fitting-encodingdecoding-models)
  - [Decoding attention and multisensory processing](#decoding-attention-and-multisensory-processing)
  - [Efficient covariance matrix estimation](#efficient-covariance-matrix-estimation)
  - [Feature extraction](#feature-extraction)
- [Examples](#examples)
  - [TRF/STRF estimation](#trfstrf-estimation)
  - [Stimulus reconstruction](#stimulus-reconstruction)
  - [Single-lag decoding analysis](#single-lag-decoding-analysis)
- [License](#license)

## Installation

Download and unzip mTRF-Toolbox to a local directory, then in the MATLAB/GNU Octave command window enter:

```matlab
addpath 'directory/mTRF-Toolbox-master'
savepath
```

## Documentation

For documentation and citation, please refer to the [mTRF-Toolbox paper](docs/Crosse_etal_FrontHumNeurosci_2016.pdf).

For usage, please see [examples](#examples) and [example M-files](examples).

## mTRF Modelling Framework

mTRF-Toolbox provides a complementary forward/backward quantitative modelling framework. A forward model, known as a temporal response function or temporal receptive field (TRF), describes how sensory information is encoded in neuronal activity. Multivariate stimulus features such as spatio- or spectro-temporal representations, as well as categorical features such as phonetic or semantic embeddings, can be used as inputs to the model. TRFs can be subjected to conventional time-frequency/source analysis techniques or used to predict the neural responses to an independent set of stimuli. mTRF-Toolbox provides an efficient cross-validation procedure for hyperparameter optimization.

A backward model, known as a neural decoder, reverses the direction of causality between stimulus and response. Neural decoders can be used to reconstruct stimulus features from information encoded explicitly or implicitly in neuronal activity, or decode higher-order cognitive processes such as top-down attention. The mTRF modelling framework provides a basic machine learning platform for real-time BCI applications such as stimulus reconstruction/synthesis and auditory attention decoding (AAD).

<div align="center">
  <img src="docs/mTRF_modelling_framework.PNG">
</div>

## Contents

### Fitting encoding/decoding models

* `mTRFcrossval()` - cross-validation for hyperparameter optimization
* `mTRFtrain()` - fits an encoding/decoding model (TRF/STRF estimation)
* `mTRFtransform()` - transforms a decoding model into an encoding model
* `mTRFpredict()` - predicts and evaluates the output of a model
* `mTRFevaluate()` - evaluates the accuracy and error of a prediction

### Decoding attention and multisensory processing

* `mTRFattncrossval()` - cross-validation for attention decoder optimization
* `mTRFattnevaluate()` - evaluates the accuracy of an attention decoder
* `mTRFmulticrossval()` - cross-validation for multisensory model optimization
* `mTRFmultitrain()` - fits an additive multisensory model (TRF/STRF estimation)

### Efficient covariance matrix estimation

* `olscovmat()` - covariance matrices for ordinary least squares estimation
* `mlscovmat()` - covariance matrices for multisensory least squares estimation

### Feature extraction

* `mTRFenvelope()` - computes the temporal envelope of an audio signal
* `mTRFresample()` - resamples and smooths temporal features
* `lagGen()` - generates time-lagged input features of multivariate data

## Examples

### TRF/STRF estimation

Here, we estimate a 16-channel spectro-temporal response function (STRF) from 2 minutes of EEG recorded while a human participant listened to natural speech.

```matlab
% Load example speech dataset
load('data/speech_data.mat','stim','resp','fs','factor');       

% Estimate STRF model weights
model = mTRFtrain(stim,resp*factor,fs,1,-100,400,0.1);
```

We can estimate the broadband TRF by taking the sum across frequency channels, and the global field power (GFP) by taking the standard deviation across EEG channels.

```
% Compute broadband TRF
strf = model.w;
trf = squeeze(sum(model.w));

% Compute global field power
sgfp = squeeze(std(strf,[],3));
gfp = std(trf,[],2);
```

We can then plot the STRF and GFP as a function of time lags. This example can also be generated using [plot_speech_STRF](examples/plot_speech_strf.m) and [plot_speech_TRF](examples/plot_speech_trf.m).

```
% Plot STRF
figure
subplot(2,2,1), imagesc(model.t(7:59),1:16,strf(:,7:59,85)), axis square
title('Speech STRF (Fz)'), ylabel('Frequency band'), set(gca,'ydir','normal')

% Plot GFP
subplot(2,2,2), imagesc(model.t(7:59),1:16,sgfp(:,7:59)), axis square
title('Global Field Power'), set(gca,'ydir','normal')

% Plot TRF
subplot(2,2,3), plot(model.t,trf(:,85),'linewidth',3), xlim([-50,350]), axis square, grid on
title('Speech TRF (Fz)'), xlabel('Time lag (ms)'), ylabel('Amplitude (a.u.)')

% Plot GFP
subplot(2,2,4), area(model.t,squeeze(gfp),'edgecolor','none'), xlim([-50,350]), axis square, grid on
title('Global Field Power'), xlabel('Time lag (ms)')
```

<img src="docs/STRF_example.PNG">

### Stimulus reconstruction

Here, we perform cross-validation (CV) to optimize the performance of a neural decoder, and then test the optimized decoder on a held-out dataset. This example can also be generated using [stimulus_reconstruction](examples/stimulus_reconstruction.m).

First, we allocate training and test sets, and then run a 10-fold CV to find the regularization value that optimizes the decoders ability to predict new stimulus features.

```matlab
% Load data
load('data/speech_data.mat','stim','resp','fs');

% Normalize and downsample data
stim = resample(sum(stim,2),64,fs);
resp = resample(resp/std(resp(:)),64,fs);
fs = 64;

% Generate training/test sets
[stimtrain,resptrain,stimtest,resptest] = mTRFpartition(stim,resp,11,1);

% Model hyperparameters
dir = -1;
tmin = 0;
tmax = 250;
lambda = 10.^(-6:2:6);
nlambda = length(lambda);

% Run fast cross-validation
cv = mTRFcrossval(stimtrain,resptrain,fs,dir,tmin,tmax,lambda,'zeropad',0,'fast',1);
```

Based on the CV analysis, we train our model using the optimal regularization value and test it on a dataset that was held-out during CV. Model performance is evaluated by measuring the correlation between the original and predicted stimulus.

``` matlab
% Use optimal regularization value
[rmax,idx] = max(mean(cv.acc));
lambda = lambda(idx);

% Train model
model = mTRFtrain(stimtrain,resptrain,fs,dir,tmin,tmax,lambda,'zeropad',0);

% Test model
[pred,test] = mTRFpredict(stimtest,resptest,model,'zeropad',0);

% Plot CV accuracy
figure
subplot(2,2,1), errorbar(1:nlambda,mean(cv.acc),std(cv.acc)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1]), axis square, grid on
title('CV Accuracy'), xlabel('Lambda (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:nlambda,mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1]), axis square, grid on
title('CV Error'), xlabel('Lambda (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(stimtest))/fs,stimtest,'linewidth',2), hold on
plot((1:length(pred))/fs,pred,'linewidth',2), hold off, xlim([0,10]), axis square, grid on
title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,rmax), hold on, bar(2,test.acc), hold off
set(gca,'xtick',1:2,'xticklabel',{'CV','Test'}), axis square, grid on
title('Test Result'), xlabel('Metric'), ylabel('Correlation')
```

<img src="docs/stim_recon_example.PNG">

### Single-lag decoding analysis

Here, we evaluate the contribution of individual time lags towards stimulus reconstruction using a single-lag decoding analysis. We perform a 10-fold cross-validation to test a series of single-lag decoders over the range 0 to 1000 ms using an optimized regularization parameter. This example can also be generated using [single_lag_analysis](examples/single_lag_analysis.m).

```matlab
% Load data
load('data/speech_data.mat','stim','resp','fs');

% Normalize and downsample data
stim = resample(sum(stim,2),64,fs);
resp = resample(resp/std(resp(:)),64,fs);
fs = 64;

% Generate training/test sets
[stimtrain,resptrain] = mTRFpartition(stim,resp,10);

% Run single-lag cross-validation
[stats,t] = mTRFcrossval(stimtrain,resptrain,fs,-1,0,1e3,10.^-2,'type','single','zeropad',0);

% Compute mean and variance
macc = squeeze(mean(stats.acc))'; vacc = squeeze(var(stats.acc))';
merr = squeeze(mean(stats.err))'; verr = squeeze(var(stats.err))';

% Compute variance bound
xacc = [-fliplr(t),-t]; yacc = [fliplr(macc-sqrt(vacc/nfold)),macc+sqrt(vacc/nfold)];
xerr = [-fliplr(t),-t]; yerr = [fliplr(merr-sqrt(verr/nfold)),merr+sqrt(verr/nfold)];

% Plot accuracy
figure
subplot(1,2,1), h = fill(xacc,yacc,'b','edgecolor','none'); hold on
set(h,'facealpha',0.2), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(t),fliplr(macc),'linewidth',2), hold off
title('Reconstruction Accuracy'), xlabel('Time lag (ms)'), ylabel('Correlation')

% Plot error
subplot(1,2,2)
h = fill(xerr,yerr,'b','edgecolor','none'); hold on
set(h,'facealpha',0.2), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(t),fliplr(merr),'linewidth',2), hold off
title('Reconstruction Error'), xlabel('Time lag (ms)'), ylabel('MSE')
```

<img src="docs/single_lag_analysis_example.PNG">

## License

[BSD 3-Clause License](LICENSE)
