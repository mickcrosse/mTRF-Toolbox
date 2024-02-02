function neural_prediction
%NEURAL_PREDICTION  Neural prediction example.
%   NEURAL_PREDICTION loads an example dataset and trains a forward
%   encoding model (TRF) that predicts 128-channel EEG responses from 2
%   minutes of speech as per Crosse et al. (2016).
%
%   Example data is loaded from SPEECH_DATA.MAT and consists of the
%   following variables:
%       'stim'      a vector containing the speech spectrogram, obtained by
%                   band-pass filtering the speech signal into 128
%                   logarithmically spaced frequency bands between 100
%                   and 4000Hz, taking the Hilbert transform at each band
%                   and averaging over every 8 neighbouring bands.
%       'resp'      a matrix containing 2 minutes of 128-channel EEG data
%                   filtered between 0.5 and 15Hz.
%       'fs'        the sample rate of STIM and RESP (128Hz).
%       'factor'    the BioSemi EEG normalization factor for converting the
%                   TRF to microvolts (524.288mV / 2^24bits).
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Crosse MJ, Zuk NJ, Di Liberto GM, Nidiffer AR, Molholm S, Lalor
%          EC (2021) Linear Modeling of Neurophysiological Responses to
%          Speech and Other Continuous Stimuli: Methodological
%          Considerations for Applied Research. Front Neurosci 15:705621.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Load data
load('../data/speech_data.mat','stim','resp','fs','factor');

% Normalize data
stim = sum(stim,2);
resp = factor*resp(:,85);

% ---Cross-validation---

% Generate training/test sets
nfold = 6;
testfold = 3;
[stimtrain,resptrain,stimtest,resptest] = mTRFpartition(stim,resp,nfold,...
    testfold);

% Model hyperparameters
Dir = -1; % direction of causality
tmin = 0; % minimum time lag
tmax = 200; % maximum time lag
lambdas = 10.^(-3:1:3); % regularization values

% Run fast cross-validation
cv = mTRFcrossval(stimtrain,resptrain,fs,Dir,tmin,tmax,lambdas,...
    'zeropad',0,'fast',1);

% ---Model training---

% Get optimal hyperparameters
[rmax,idx] = max(mean(mean(cv.r),3));
lambda = lambdas(idx);
nlambda = length(lambdas);

% Train model
model = mTRFtrain(stimtrain,resptrain,fs,Dir,tmin,tmax,lambda,'zeropad',0);

% ---Model testing---

% Test model
[pred,test] = mTRFpredict(stimtest,resptest,model,'zeropad',0);

% ---Evaluation---

% Plot CV correlation
figure('Name','Neural Prediction','NumberTitle','off')
set(gcf,'color','w')
subplot(2,2,1)
errorbar(1:nlambda,mean(cv.r),std(cv.r)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1])
title('CV Accuracy')
xlabel('Regularization (1\times10^\lambda)')
ylabel('Correlation')
axis square, grid on

% Plot CV error
subplot(2,2,2)
errorbar(1:nlambda,mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1])
title('CV Error')
xlabel('Regularization (1\times10^\lambda)')
ylabel('MSE')
axis square, grid on

% Plot reconstruction
subplot(2,2,3)
plot((1:length(resptest))/fs,resptest,'linewidth',2), hold on
plot((1:length(pred))/fs,pred,'linewidth',2), hold off
xlim([0,10])
title('Prediction')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
axis square, grid on
legend('Orig','Pred')

% Plot test correlation
subplot(2,2,4)
bar(1,rmax), hold on
bar(2,test.r), hold off
xlim([0,3])
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'})
title('Model Performance')
xlabel('Dataset')
ylabel('Correlation')
axis square, grid on