function plot_contrast_trf
%PLOT_CONTRAST_TRF  Plot example contrast TRF (VESPA).
%   PLOT_CONTRAST_TRF loads an example dataset, estimates and plots a
%   contrast TRF (VESPA) and the global field power (GFP) from 2 minutes of
%   128-channel EEG data as per Lalor et al. (2006).
%
%   Example data is loaded from CONTRAST_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the normalized contrast levels of a
%                   checkerboard that modulated at a rate of 60Hz
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 35 Hertz
%       'fs'        the sample rate of STIM and RESP (128Hz)
%       'Nf'        the modulation rate of the checkerboard (60Hz)
%       'factor'    the BioSemi EEG normalization factor for computing the
%                   TRF in microvolts (524.288mV / 2^24bits)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Lalor EC, Pearlmutter BA, Reilly RB, McDarby G, Foxe JJ (2006)
%          The VESPA: a method for the rapid estimation of a visual evoked
%          potential. NeuroImage 32:1549-1561.
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
load('../data/contrast_data.mat','stim','resp','fs','Nf','factor');

% Normalize data
stim = stim*Nf;
resp = resp*factor;

% Model hyperparameters
Dir = 1; % direction of causality
tmin = -100; % minimum time lag
tmax = 400; % maximum time lag
lambda = 0.5; % regularization value

% Compute model weights
model = mTRFtrain(stim,resp,fs,Dir,tmin,tmax,lambda,'method','Tikhonov',...
    'split',5,'zeropad',0);

% Plot TRF
figure('Name','Contrast TRF & GFP','NumberTitle','off')
set(gcf,'color','w')
subplot(1,2,1)
mTRFplot(model,'trf',[],23,[-50,350]);
title('Contrast TRF (Oz)')
ylabel('Amplitude (\muV)')

% Plot GFP
subplot(1,2,2)
mTRFplot(model,'gfp',[],'all',[-50,350]);
title('Global Field Power')