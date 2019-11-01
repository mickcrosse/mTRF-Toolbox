% Demo for speech data

addpath('..');

lambdas = 10.^(-2:0.5:4);
minlag = -150; % minimum lag in the model (ms)
maxlag = 450; % maximum lag in the model (ms)

% Load the speech data
load('speech_data.mat');

%% Envelope model
disp('** Envelope **');
% Train the model using 10-fold cross-validation
[r_env,~,~,env_forward_model] = mTRFcrossval(envelope,EEG,Fs,1,minlag,maxlag,lambdas);

% Identify the optimal lambda value
R_env = mean(r_env,3); % average r across channels
opt_lmb_idx = find(mean(R_env)==max(mean(R_env)),1);
fprintf('Optimal lambda = %.2f\n',lambdas(opt_lmb_idx));
opt_env_model = env_forward_model(:,:,:,opt_lmb_idx);

% Plot the TRF model
nchans = size(EEG,2);
plot_trf('Speech envelope',opt_env_model,Fs,minlag,maxlag);

%% Spectrogram model
disp('** Spectrogram **');
% Train the model using 10-fold cross-validation
[r_spec,~,~,spec_forward_model] = mTRFcrossval(spectrogram,EEG,Fs,1,minlag,maxlag,lambdas);

% Identify the optimal lambda value
R_spec = mean(r_spec,3); % average r across channels
opt_lmb_idx = find(mean(R_spec)==max(mean(R_spec)),1);
fprintf('Optimal lambda = %.2f\n',lambdas(opt_lmb_idx));
opt_spec_model = spec_forward_model(:,:,:,opt_lmb_idx);

% Plot the mTRF spectrogram model, 16 different frequencies
plot_multifeature_trf('Speech spectrogram',opt_spec_model,Fs,minlag,maxlag,1:16);