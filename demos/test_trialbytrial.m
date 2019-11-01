% Do trail-by-trial testing for a forward model

addpath('..');

ntrials = 4; % # of trials in dataset
sbj = 'Subject10'; % subject label
fs = 128; % sampling rate (in Hz) of the stimulus envelope and EEG in the dataset
data_path = '..\sample_data\';
trf_maxlag = 350; % maximum lag in the TRF (in ms)
trf_minlag = -50; % minimum lag in the TRF (in ms)
% lambda = 10.^(0:0.5:3); % use a single lambda value
lambda = 10.^5;

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

%% Do leave-one-out testing
[r,p,rmse,yhat,w,t,b,opt_lambda] = mTRFtrialbytrial(stims,eegs,fs,-1,trf_minlag,trf_maxlag,lambda);

%% Plotting
% Plot one of the models