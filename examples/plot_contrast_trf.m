function plot_contrast_trf
%PLOT_CONTRAST_TRF  Plot example contrast TRF (VESPA).
%   PLOT_CONTRAST_TRF loads an example dataset and estimates and plots a
%   contrast TRF (VESPA) and the global field power (GFP) from 2 minutes of
%   128-channel EEG data as per Lalor et al. (2006).
%
%   Example data is loaded from CONTRAST_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the normalized contrast levels of a
%                   checkerboard that modulated at a rate of 60Hz
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 35 Hz
%       'fs'        the sample rate of STIM and RESP (128Hz)
%       'Nf'        the modulation rate of the checkerboard (60Hz)
%       'factor'    the BioSemi EEG normalization factor for converting the
%                   TRF to microvolts (524.288mV / 2^24bits)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Lalor EC, Pearlmutter BA, Reilly RB, McDarby G, Foxe JJ (2006)
%          The VESPA: a method for the rapid estimation of a visual evoked
%          potential. NeuroImage 32:1549-1561.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Load data
load('data/contrast_data.mat','stim','resp','fs','Nf','factor');

% Normalize data to convert TRF to uV
stim = stim*Nf;
resp = resp*factor;

% Model hyperparameters
dir = 1;
tmin = -100;
tmax = 400;
lambda = 4.4e-3;

% Compute model weights
model = mTRFtrain(stim,resp,fs,dir,tmin,tmax,lambda,'method','Tikhonov',...
    'zeropad',0);

% Get TRF and GFP
trf = squeeze(model.w);
gfp = squeeze(std(model.w,[],3));

% Define ROI
chan = 23; % channel Oz

% Plot TRF
figure, subplot(1,2,1)
plot(model.t,trf(:,chan),'linewidth',3)
xlim([-50,350])
title('Contrast TRF (Oz)')
xlabel('Time lag (ms)')
ylabel('Amplitude (\muV)')
axis square, grid on

% Plot GFP
subplot(1,2,2)
area(model.t,gfp,'edgecolor','none');
xlim([-50,350])
title('Global Field Power')
xlabel('Time lag (ms)')
axis square, grid on