% EEG prediction demo, based on envelope, show tuning curves as a function of lambda
% using 10-fold cross-validation and trial-by-trial cross-validation
% The data for this demo is taken from the dryad dataset by Broderick et
% al (2018): https://doi.org/10.1016/j.cub.2018.01.080
% Nate Zuk (2019)

addpath('..');

ntrials = 10; % # of trials in dataset
sbj = 'Subject2'; % subject label
fs = 128; % sampling rate (in Hz) of the stimulus envelope and EEG in the dataset
data_path = '..\sample_data\';
trf_maxlag = 350; % maximum lag in the TRF (in ms)
trf_minlag = -50; % minimum lag in the TRF (in ms)
lambdas = 10.^(0:0.5:10);

% Load the stimulus envelopes
disp('Loading the stimulus envelopes...');
stims = cell(ntrials,1); % setup the cell array to store envelopes
for n = 1:ntrials
    % get the filename for the stimulus file
    stim_filename = sprintf('audio%d_128Hz',n);
    % load the stimulus envelope file
    S = load([data_path 'envs\' stim_filename]);
    % store the envelope
    stims{n} = S.env;
    disp(stim_filename); % display filename to show successful loading
end
clear S

% Load the EEG data
disp('Loading the EEG data...');
eegs = cell(ntrials,1); % set up the cell array to store EEG
for n = 1:ntrials
    % get the filename for the stimulus file
    eeg_filename = sprintf('%s_Run%d',sbj,n);
    % load the stimulus envelope file
    S = load([data_path 'eegs\' eeg_filename]);
    % remove the average of the mastoids, and store the EEG
    eegs{n} = S.eegData-mean(S.mastoids,2)*ones(1,128);
    disp(eeg_filename); % display filename to show successful loading
end
clear S

%%% Do very basic EEG preprocessing -- highpass filter the EEG using a zero-phase 
%%% moving average filter with a window equal to the size of the TRF
disp('Do some very basic EEG preprocessing (high pass filter...');
% compute the window size
trf_window_size = ceil(trf_maxlag/1000*fs)-floor(trf_minlag/1000*fs);
fprintf('...using a moving average window of length %.0f ms)\n',trf_window_size/fs*1000);
for n = 1:ntrials
    % remove the moving average with that window size
    eegs{n} = eegs{n}-movmean(eegs{n},trf_window_size);
end

fprintf('\n'); % put a break in the output text, before the modeling starts

%%% Run the mTRFcrossval with map=-1 (backwards modeling, envelope reconstruction)
disp('** Doing 10-fold cross-validation, with a random sampling of data for each fold **');
[r_cv,~,~,model_cv,t] = mTRFcrossval(stims,eegs,fs,1,trf_minlag,trf_maxlag,lambdas);
% average across channels
r_cv = mean(r_cv,3);

%%% Run leave-one-out, trial-by-trial cross-validation
disp('** Doing leave-one-out, trial-by-trial cross-validation **');
[r_loo,~,~,model_loo] = mTRFcrossval_loo(stims,eegs,fs,1,trf_minlag,trf_maxlag,lambdas);
r_loo = mean(r_loo,3);

%%% Plot the reconstruction accuracy as a function of lambda, for both r
%%% and mse
figure
% plot r for 10-fold cross-validation
subplot(1,2,1);
plot(lambdas,r_cv);
set(gca,'XScale','log','FontSize',16);
xlabel('\lambda');
ylabel('Pearson''s correlation');
opt_lambda_rcv = lambdas(mean(r_cv)==max(mean(r_cv))); % compute the optimum lambda
tle = sprintf('10-fold CV (opt lambda = %.0f)',opt_lambda_rcv);
title(tle);

% plot r for leave-one-out cross-validation
subplot(1,2,2);
plot(lambdas,r_loo);
set(gca,'XScale','log','FontSize',16);
xlabel('\lambda');
ylabel('Pearson''s correlation');
opt_lambda_rloo = lambdas(mean(r_loo)==max(mean(r_loo))); % compute the optimum lambda
tle = sprintf('Leave-one-out CV (opt lambda = %.0f)',opt_lambda_rloo);
title(tle);