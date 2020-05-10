function plot_speech_strf
%PLOT_SPEECH_STRF  Plot example speech STRF.
%   PLOT_SPEECH_STRF loads an example dataset, estimates and plots a speech
%   STRF and the global field power (GFP) from 2 minutes of 128-channel EEG
%   data as per Di Liberto et al. (2015).
%
%   Example data is loaded from SPEECH_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 15 Hertz
%       'fs'        the sample rate of STIM and RESP (128Hz)
%       'factor'    the BioSemi EEG normalization factor for computing the
%                   STRF in microvolts (524.288mV / 2^24bits)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Di Liberto GM, O'Sullivan JA, Lalor EC (2015) Low-Frequency
%          Cortical Entrainment to Speech Reflects Phoneme-Level
%          Processing. Curr Biol 25:1-9.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Load data
load('data/speech_data.mat','stim','resp','fs','factor');

% Normalize data
resp = factor*resp;

% Model hyperparameters
tmin = -100;
tmax = 400;
lambda = 0.5;

% Compute model weights
% Note, ridge regression is used instead of Tikhonov regularization to
% avoid cross-channel leakage of the multivariate input features
model = mTRFtrain(stim,resp,fs,1,tmin,tmax,lambda,'method','ridge',...
    'split',5,'zeropad',0);

% Plot STRF
figure, subplot(1,2,1)
mTRFplot(model,'mtrf','all',85,[-50,350]);
title('Speech STRF (Fz)')
ylabel('Frequency band')

% Plot GFP
subplot(1,2,2)
mTRFplot(model,'mgfp','all','all',[-50,350]);
title('Global Field Power')