function stimulus_reconstruction
%STIMULUS_RECONSTRUCTION  Stimulus reconstruction example.
%   STIMULUS_RECONSTRUCTION loads an example dataset and trains a neural
%   decoder that reconstructs stimulus features (speech envelope) from 2
%   minutes of 128-channel EEG data as per Crosse et al. (2016).
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
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Load data
load('data/speech_data.mat','stim','resp','fs');

% Normalize data
stim = sum(stim,2);
resp = resp/std(resp(:));

% Downsample data
fsNew = 64;
stim = resample(stim,fsNew,fs);
resp = resample(resp,fsNew,fs);
fs = fsNew;

% ---Create training/test sets---

%TODO(mickcrosse): create cvfold() function

% Define cross-validation folds
nfold = 10;
batch = ceil(size(stim,1)/nfold);

% Training set
stimtrain = cell(nfold,1);
resptrain = cell(nfold,1);
for i = 1:nfold
    idx = batch*(i-1)+1:min(batch*i,size(resp,1));
    stimtrain{i} = stim(idx,:);
    resptrain{i} = resp(idx,:);
end

% Test set
itest = 1;
stimtest = stimtrain{itest};
resptest = resptrain{itest};

% Remove test set from training data
stimtrain(itest) = [];
resptrain(itest) = [];

% ---Cross-validation---

% Model hyperparameters
dir = -1;
tmin = 0;
tmax = 250;
lambdas = 10.^(-6:2:6);
nlambda = length(lambdas);

% Run fast cross-validation
cv = mTRFcrossval(stimtrain,resptrain,fs,dir,tmin,tmax,lambdas,...
    'zeropad',0,'fast',1);

% ---Model training---

% Get optimal hyperparameters
[rmax,idx] = max(mean(cv.acc));
lambda = lambdas(idx);

% Train model
model = mTRFtrain(stimtrain,resptrain,fs,dir,tmin,tmax,lambda,...
    'zeropad',0);

% ---Model testing---

% Test model
[pred,test] = mTRFpredict(stimtest,resptest,model,'zeropad',0);

% ---Evaluation---

% Plot CV accuracy
figure
subplot(2,2,1)
errorbar(1:nlambda,mean(cv.acc),std(cv.acc)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1])
title('CV Accuracy')
xlabel('Lambda (1\times10^\lambda)')
ylabel('Correlation')
axis square, grid on

% Plot CV error
subplot(2,2,2)
errorbar(1:nlambda,mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6), xlim([0,nlambda+1])
title('CV Error')
xlabel('Lambda (1\times10^\lambda)')
ylabel('MSE')
axis square, grid on

% Plot reconstruction
subplot(2,2,3)
plot((1:length(stimtest))/fs,stimtest,'linewidth',2), hold on
plot((1:length(pred))/fs,pred,'linewidth',2), hold off
xlim([0,10])
title('Reconstruction')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
axis square, grid on
legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4)
bar(1,rmax), hold on
bar(2,test.acc), hold off
set(gca,'xtick',1:2,'xticklabel',{'CV','Test'})
title('Test Result')
xlabel('Metric')
ylabel('Correlation')
axis square, grid on