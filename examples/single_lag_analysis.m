function single_lag_analysis
%SINGLE_LAG_ANALYSIS  Single-lag analysis example.
%   SINGLE_LAG_ANALYSIS loads an example dataset and trains a set of
%   single-lag neural decoders that reconstruct stimulus features (speech
%   envelope) from 2 minutes of 128-channel EEG data as per O'Sullivan et
%   al (2015).
%
%   Example data is loaded from SPEECH_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 15 Hz
%       'fs'        the sample rate of STIM and RESP (128Hz)
%       'factor'    the BioSemi EEG normalization factor for converting the
%                   TRF to microvolts (524.288mV / 2^24bits)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] O'Sullivan JA, Power AJ, Mesgarani N, Rajaram S, Foxe JJ, Shinn-
%          Cunningham BG, Slaney M, Shamma SA, Lalor EC (2015) Attentional
%          Selection in a Cocktail Party Environment Can Be Decoded from
%          Single-Trial EEG. Cereb Cortex 25(7):1697-1706.
%      [2] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [3] Crosse MJ, Zuk NJ, Di Liberto GM, Nidiffer AR, Molholm S, Lalor 
%          EC (2021) Linear Modeling of Neurophysiological Responses to 
%          Speech and Other Continuous Stimuli: Methodological 
%          Considerations for Applied Research. Front Neurosci 15:705621.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Load data
load('../data/speech_data.mat','stim','resp','fs');

% Normalize data
stim = sum(stim,2);
resp = resp/std(resp(:));

% Downsample data
fsNew = 64;
stim = resample(stim,fsNew,fs);
resp = resample(resp,fsNew,fs);
fs = fsNew;

% Generate training/test sets
nfold = 10;
[stimtrain,resptrain] = mTRFpartition(stim,resp,nfold);

% Model hyperparameters
Dir = -1; % direction of causality
tmin = 0; % minimum time lag
tmax = 1000; % maximum time lag
lambda = 10.^-2; % regularization value

% Run single-lag cross-validation
[stats,t] = mTRFcrossval(stimtrain,resptrain,fs,Dir,tmin,tmax,...
    lambda,'type','single','zeropad',0);

% Compute mean and variance
mr = squeeze(mean(stats.r))';
merr = squeeze(mean(stats.err))';
vr = squeeze(var(stats.r))';
verr = squeeze(var(stats.err))';

% Compute variance bound
xr = [-fliplr(t),-t];
xerr = [-fliplr(t),-t];
yr = [fliplr(mr-sqrt(vr/nfold)),mr+sqrt(vr/nfold)];
yerr = [fliplr(merr-sqrt(verr/nfold)),merr+sqrt(verr/nfold)];

% Plot accuracy
figure('Name','Single-lag Analysis','NumberTitle','off')
set(gcf,'color','w')
subplot(1,2,1)
h = fill(xr,yr,'b','edgecolor','none'); hold on
set(h,'facealpha',0.2), xlim([tmin,tmax])
plot(-fliplr(t),fliplr(mr),'linewidth',2), hold off
title('Reconstruction Accuracy')
xlabel('Time lag (ms)')
ylabel('Correlation')
axis square, grid on

% Plot error
subplot(1,2,2)
h = fill(xerr,yerr,'b','edgecolor','none'); hold on
set(h,'facealpha',0.2), xlim([tmin,tmax])
plot(-fliplr(t),fliplr(merr),'linewidth',2), hold off
title('Reconstruction Error')
xlabel('Time lag (ms)')
ylabel('MSE')
axis square, grid on
