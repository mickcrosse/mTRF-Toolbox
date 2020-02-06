% mTRF-Toolbox demos.
%
% This script contains several examples to demonstrate the functionality of
% mTRF-Toolbox. To run each example individually, place the cursor in the 
% seciton of interest and click the 'Run Section' button in the EDITOR tab,
% or press Ctrl + Enter (Windows) or Cmd + Return (Mac).

% URL: https://github.com/mickcrosse/mTRF-Toolbox

%% 1. Visual Contrast: TRF Estimation

clear;

% Load data
load('contrast_data.mat');

% Normalize data
stim = contrastLevel/std(contrastLevel); 
resp = EEG/std(EEG(:));

% Model parameters
dir = 1;
tmin = -150;
tmax = 450;
lambda = 0;

% Compute model weights
model = mTRFtrain(stim,resp,Fs,dir,tmin,tmax,lambda,'zeropad',0);

% Define ROI 
chan = 23; % channel Oz

% Plot TRF & GFP
figure(1)
subplot(1,2,1)
plot(model.t,squeeze(model.w(:,:,chan)),'linewidth',3)
xlim([-50,350])
title ('Contrast TRF (Oz)')
xlabel('Time lag (ms)')
ylabel('Amplitude (a.u.)')
axis square
grid on
subplot(1,2,2)
area(model.t,squeeze(std(model.w,[],3)),'edgecolor','none');
xlim([-50,350])
title ('Global Field Power')
xlabel('Time lag (ms)')
axis square
grid on

%% 2. Audio Speech: TRF Estimation

clear;

% Load data
load('speech_data.mat');

% Normalize data
stim = sum(spectrogram,2); 
stim = stim/std(stim); 
resp = EEG/std(EEG(:));

% Model parameters
dir = 1;
tmin = -150;
tmax = 450;
lambda = 10;

% Compute model weights
model = mTRFtrain(stim,resp,Fs,dir,tmin,tmax,lambda,'method','tik',...
    'zeropad',0);
    
% Define ROI 
chan = 85; % channel Fz

% Plot TRF & GFP
figure(2)
subplot(1,2,1)
plot(model.t,squeeze(model.w(:,:,chan)),'linewidth',3)
xlim([-50,350]);
title ('Speech TRF (Fz)')
xlabel('Time lag (ms)')
ylabel('Amplitude (a.u.)')
axis square
grid on
subplot(1,2,2)
area(model.t,squeeze(std(model.w,[],3)),'edgecolor','none');
xlim([-50,350])
title ('Global Field Power')
xlabel('Time lag (ms)')
axis square
grid on

%% 3. Audio Speech: STRF Estimation

clear; clc;

% Load data
load('speech_data.mat');

% Normalize data
stim = spectrogram/std(spectrogram(:)); 
resp = EEG/std(EEG(:));
n = size(stim,2);

% Model parameters
dir = 1;
tmin = -150;
tmax = 450;
lambda = 10*n^2;

% Compute model weights
model = mTRFtrain(stim,resp,Fs,dir,tmin,tmax,lambda,'zeropad',0);

% Define ROI 
chan = 85; % channel Fz

% Plot STRF & GFP
figure(3)
subplot(1,2,1)
imagesc(model.t(14:66),1:16,squeeze(model.w(:,14:66,chan)))
xlim([-50,350])
set(gca,'ydir','normal')
title ('Speech STRF (Fz)')
xlabel('Time lag (ms)')
ylabel('Frequency band')
axis square
subplot(1,2,2)
imagesc(model.t(14:66),1:16,squeeze(std(model.w(:,14:66,:),[],3)))
xlim([-50,350])
set(gca,'ydir','normal')
title ('Global Field Power')
xlabel('Time lag (ms)')
axis square

% Compute broadband TRF by summing over STRF channels:

% TRF parameters
dir = 1;
tmin = -150;
tmax = 450;
lambda = 10*n;

% Compute model weights
model = mTRFtrain(stim,resp,Fs,dir,tmin,tmax,lambda,'method','rid');

% Define ROI 
chan = 85; % channel Fz

% Plot TRF & GFP
figure(4)
subplot(1,2,1)
plot(model.t,squeeze(sum(model.w(:,:,chan))),'linewidth',3)
xlim([-50,350])
title ('Broadband TRF (Fz)')
xlabel('Time lag (ms)')
ylabel('Amplitude (a.u.)')
axis square
grid on
subplot(1,2,2)
area(model.t,squeeze(std(sum(model.w),[],3)),'edgecolor','none')
xlim([-50,350])
title ('Global Field Power')
xlabel('Time lag (ms)')
axis square
grid on

%% 4. Audio Speech: Stimulus Reconstruction

clear; clc;

% Load data
load('speech_data.mat');

% Normalize data
stim = sum(spectrogram,2); 
stim = stim/std(stim); 
resp = EEG/std(EEG(:));

% Downsample data
FsNew = 64;
stim = resample(stim,FsNew,Fs); 
resp = resample(resp,FsNew,Fs);
Fs = FsNew;

% Determine folds for cross-validation
nFolds = 10;
batch_size = floor(size(stim,1)/nFolds);

% Training set
k = 1;
stimtrain = cell(nFolds,1);
resptrain = cell(nFolds,1);
for i = 1:nFolds
    idx = k:k+batch_size-1;
    stimtrain{i} = stim(idx,:); 
    resptrain{i} = resp(idx,:);
    k = k+batch_size;
end
    
% Test set
% itest = randi([1,nFolds],1); % pick test set at random
itest = 1;
stimtest = stimtrain{itest}; 
resptest = resptrain{itest}; 
stimtrain(itest) = []; % remove test set from training data
resptrain(itest) = []; 

% % Reduce EEG dimensions using PCA
% npcs = 64;
% [resptrain,V] = mTRFpca(resptrain,npcs,'pca','eig');
% resptest = resptest*V;

% Cross-validation parameters
dir = -1;
tmin = 50;
tmax = 200;
lambdas = 10.^(-6:2:6);
nlambda = length(lambdas);

% Run cross-validation
[r,~,rmse] = mTRFcrossval(stimtrain,resptrain,Fs,dir,tmin,tmax,lambdas,...
    'zeropad',0,'fast',1);        

% Plot cross-validation results
figure(5)
subplot(1,2,1)
errorbar(1:nlambda,mean(r),std(r)/sqrt(nFolds-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6)
xlim([0,nlambda+1])
title ('Cross-Validation (Acc.)')
xlabel('Regularization (1\times10^\lambda)')
ylabel('Correlation')
axis square
grid on
subplot(1,2,2)
errorbar(1:nlambda,mean(rmse),std(rmse)/sqrt(nFolds-1),'linewidth',2)
set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6)
xlim([0,nlambda+1])
title ('Cross-Validation (Err.)')
xlabel('Regularization (1\times10^\lambda)')
ylabel('RMSE')
axis square
grid on

% Model parameters
[rmax,idx] = max(mean(r));
lambda = lambdas(idx);

% Train model
model = mTRFtrain(stimtrain,resptrain,Fs,dir,tmin,tmax,lambda,'zeropad',1);

% Test model
[pred,stats] = mTRFpredict(stimtest,resptest,model,'zeropad',1);

% Plot test results
figure(6)
subplot(1,2,1)
plot(stimtest,'linewidth',2)
hold on
plot(pred,'linewidth',2)
hold off
xlim([0,400])
title ('Stimulus Reconstruction')
xlabel('Sample number')
ylabel('Amplitude (a.u.)')
axis square
grid on
legend('Orig','Pred')
subplot(1,2,2)
bar(1,rmax)
hold on
bar(2,stats.r)
title ('Test Result')
xlabel('Metric')
ylabel('Correlation')
set(gca,'xtick',1:2,'xticklabel',{'CV mean','Test'})
hold off
axis square
grid on

%% 5. Audio Speech: Transform backward model to forward model

clear; clc;

% Load data
load('speech_data.mat');

% Normalize data
stim = sum(spectrogram,2); 
stim = stim/std(stim); 
resp = EEG/std(EEG(:)); 

% Model parameters
dir = -1;
tmin = -150;
tmax = 450;
lambda = 50;

% Compute backward model weights
bmodel = mTRFtrain(stim,resp,Fs,dir,tmin,tmax,lambda);

% Transform to forward model weights
fmodel = mTRFtransform(bmodel,resp);

% Define ROI 
chan = 85; % channel Fz

% Plot TRF & GFP
figure(7)
subplot(1,2,1)
plot(fmodel.t,squeeze(fmodel.w(chan,:)),'linewidth',3); 
xlim([-50,350]);
title ('Frontal TRF (Fz)')
xlabel('Time lag (ms)')
ylabel('Amplitude (a.u.)')
axis square
grid on
subplot(1,2,2)
area(fmodel.t,squeeze(std(fmodel.w)),'edgecolor','none');
xlim([-50,350])
title ('Global Field Power')
xlabel('Time lag (ms)')
axis square
grid on

%% 7. Audio Speech: Sinigle-Lag Stimulus Reconstruction

clear; clc;

% Load data
load('speech_data.mat');

% Normalize data
stim = sum(spectrogram,2); 
stim = stim/std(stim(:)); 
resp = EEG/std(EEG(:));

% Downsample data
FsNew = 64;
stim = resample(stim,FsNew,Fs); 
resp = resample(resp,FsNew,Fs);
Fs = FsNew;
nobs = size(stim,1);

% Reduce EEG dimensions using PCA
npcs = 64;
[resp,scree] = mTRFpca(resp,npcs,'algorithm','eig');

% Determine folds for cross-validation
nFolds = 5;
batch_size = floor(nobs/nFolds);

% Training set
k = 1;
stimtrain = cell(nFolds,1);
resptrain = cell(nFolds,1);
for i = 1:nFolds
    idx = k:k+batch_size-1;
    stimtrain{i} = stim(idx,:); 
    resptrain{i} = resp(idx,:);
    k = k+batch_size;
end
     
% % Test set
% % itest = randi([1,nFolds],1); % pick test set at random
% itest = nFolds; % pick test set at random
% stimtest = stimtrain{itest}; 
% resptest = resptrain{itest}; 
% stimtrain(itest) = []; % remove test set from training data
% resptrain(itest) = []; 

% Cross-validation parameters
dir = -1;
tmin = -500;
tmax = 1500;
lambdas = 10.^(-6:2:6);
nlambda = length(lambdas);

% Run cross-validation
[r,~,rmse,t] = mTRFcrossval(stimtrain,resptrain,Fs,dir,tmin,tmax,lambdas,'type','sing');

% Plot reconstruction accuracy
mr = squeeze(mean(r));
vr = squeeze(var(r));
figure(8)
subplot(1,2,1)
% error = [fliplr(mr-sqrt(vr/nFolds-1)),mr+sqrt(vr/nFolds-1)];
% x = [-fliplr(t),-t];
% h = fill(x,error,'b','edgecolor','none');
% set(h,'facealpha',0.2)
% hold on
plot(-fliplr(t),fliplr(mr),'linewidth',2)
hold off
xlim([tmin,tmax])
title ('Reconstruction Accuracy')
xlabel('Time lag (ms)')
ylabel('Correlation')
axis square
grid on
subplot(1,2,2)
mrmse = mean(squeeze(rmse));
vrmse = var(squeeze(rmse));
error = [fliplr(mrmse-sqrt(vrmse/nFolds-1)),mrmse+sqrt(vrmse/nFolds-1)];
x = [-fliplr(t),-t];
h = fill(x,error,'b','edgecolor','none');
set(h,'facealpha',0.2)
hold on
plot(-fliplr(t),fliplr(mrmse),'linewidth',2)
hold off
xlim([tmin,tmax])
title ('Reconstruction Error')
xlabel('Time lag (ms)')
ylabel('RMSE')
axis square
grid on

%% 8. Downsample data

clear; clc;

% Load data
load('speech_data.mat');

% Normalize data
stim1 = sum(spectrogram,2); 
t1 = (0:length(stim1)-1)/Fs;

% Downsample data by non-integer value
Fs2 = 50;
[stim2,t2] = mTRFresample(stim1,Fs,Fs2);

% Upsample data by non-integer value
Fs3 = 200;
[stim3,t3] = mTRFresample(stim1,Fs,Fs3);

% Smooth data across 10 output frames
window = 10;
stim4 = mTRFresample(stim1,Fs,[],window);

% Plot signals
figure(9)
subplot(3,1,1)
plot(t1,stim1)
hold on
plot(t2,stim2)
hold off
subplot(3,1,2)
plot(t1,stim1)
hold on
plot(t3,stim3)
hold off
subplot(3,1,3)
plot(t1,stim1)
hold on
plot(t1,stim4)
hold off
