% Demonstrate using tlims to skip silences when computing a model for
% envelope reconstruction
% The data for this demo is taken from the dryad dataset by Broderick et
% al (2018): https://doi.org/10.1016/j.cub.2018.01.080
% !! This demo can take 20-25 minutes to complete
% Nate Zuk (2019)

addpath('..');

train_trials = 1:4; % trials to use for training
test_trials = 5; % trial to use for testing
sbj = 'Subject2'; % subject identifier
fs = 128; % sampling rate (in Hz) of the stimulus envelope and EEG in the dataset
data_path = '..\sample_data\';
trf_maxlag = 350; % maximum lag in the TRF (in ms)
trf_minlag = -50; % minimum lag in the TRF (in ms)
lambdas = 10.^(0:0.5:10);
db_thres = -67; % threshold in dB V for identifying a silent section

% Get the set of all trials so that they can be loaded in one go
all_trials = unique([train_trials test_trials]); % get the set of all training and testing trials
[~,train_idx] = intersect(all_trials,train_trials); % get the indexes in "all_trials" for training
[~,test_idx] = intersect(all_trials,test_trials); % get the indexes in "all_trials" for testing
ntrials = length(all_trials);

%% Load data
% Load the stimulus envelopes
disp('Loading the stimulus envelope...');
stims = cell(ntrials,1);
for n = 1:ntrials
    % get the filename for the stimulus file
    stim_filename = sprintf('audio%d_128Hz',all_trials(n));
    % load the stimulus envelope file
    S = load([data_path 'envs\' stim_filename]);
    % store the envelope
    stims{n} = S.env;
    disp(stim_filename); % display filename to show successful loading
end

% Load the EEG data
disp('Loading the EEG data...');
eegs = cell(ntrials,1);
for n = 1:ntrials
    % get the filename for the stimulus file
    eeg_filename = sprintf('%s_Run%d',sbj,all_trials(n));
    % load the stimulus envelope file
    S = load([data_path 'eegs\' eeg_filename]);
    % remove the average of the mastoids, and store the EEG
    eegs{n} = S.eegData-mean(S.mastoids,2)*ones(1,128);
    disp(eeg_filename); % display filename to show successful loading
end

%% EEG preprocessing
% Do very basic EEG preprocessing -- highpass filter the EEG using a zero-phase 
% moving average filter with a window equal to the size of the TRF
disp('Do some very basic EEG preprocessing (high pass filter...');
% compute the window size
trf_window_size = ceil(trf_maxlag/1000*fs)-floor(trf_minlag/1000*fs);
fprintf('...using a moving average window of length %.0f ms)\n',trf_window_size/fs*1000);
for n = 1:ntrials
    eegs{n} = eegs{n}-movmean(eegs{n},trf_window_size);
end

%% Identify silences
% Identify silent periods, and select the ranges of time
% around them
% make sure the envelope is at least 0
disp('Identifying silent periods in the stimulus...');
% concatenate stimuli together
all_stim = [];
for n = 1:ntrials, all_stim = [all_stim; stims{n}]; end
all_stim(all_stim<0) = 0;
% Get a histogram of all of the envelope values in dB, in order to identify
% the threshold for getting silences.  Usually there is quite a bit of time
% when silences occur during speech, so this distribution will be
% noticeably bimodal.
db_stim = 20*log10(all_stim); % convert to dB V
db_bins = -120:1:0;
db_histogram = histcounts(db_stim,db_bins);
clear all_stim

% plot the histogram of the envelope values in dB, and mark the threshold
% subplot(1,2,2);
figure
hold on
bar(db_bins(1:end-1)+diff(db_bins)/2,db_histogram,1,'k');
plot([db_thres db_thres],[0 max(db_histogram)],'r--','LineWidth',2);
set(gca,'FontSize',16);
xlabel('dB V');
ylabel('# samples');
legend('dB envelope values','Threshold for silences');

% Get start and end times encompassing silent periods
sil_chk = cell(ntrials,1);
sil_tlims = cell(ntrials,1);
for n = 1:ntrials
    sil_chk{n} = stims{n}<10.^(db_thres/20); % identify silences
    sil_starts = find(diff(sil_chk{n})==1); % get the indexes just before silences
    sil_ends = find(diff(-sil_chk{n})==1)+1; % get the indexes just after silences
    sil_idx = sort([sil_starts; sil_ends]); % put these indexes in ascending order
    sil_tlims{n} = (sil_idx-1)/fs; % convert to time (usetinds assumes first index is t=0)
end

fprintf('\n'); % add a break in the text before modeling

%% Training
% Do envelope reconstruction using the full envelope
disp('** Compute the model for the full envelope **');
[r_full,~,~,model_full,dly] = mTRFcrossval(stims(train_idx),eegs(train_idx),fs,-1,...
    trf_minlag,trf_maxlag,lambdas);
% get the optimal lambda
opt_lmb_idx_full = find(mean(r_full)==max(mean(r_full)),1);
fprintf('Optimal lambda = %.0f\n',lambdas(opt_lmb_idx_full));

%%% Do envelope reconstruction without silent periods
disp('** Compute the model for the envelope without silences **');
[r_nosil,~,~,model_nosil] = mTRFcrossval(stims(train_idx),eegs(train_idx),fs,-1,...
    trf_minlag,trf_maxlag,lambdas,sil_tlims(train_idx));
% get the optimal lambda
opt_lmb_idx_nosil = find(mean(r_nosil)==max(mean(r_nosil)),1);
fprintf('Optimal lambda = %.0f\n',lambdas(opt_lmb_idx_full));

%% Testing
disp('** Testing **');
% Generate the predictions for both models
[recon_full,r_full] = mTRFpredict(stims(test_idx),eegs(test_idx),...
    model_full(:,:,opt_lmb_idx_full),fs,-1,trf_minlag,trf_maxlag);
[recon_nosil,r_nosil] = mTRFpredict(stims(test_idx),eegs(test_idx),...
    model_nosil(:,:,opt_lmb_idx_nosil),fs,-1,trf_minlag,trf_maxlag,...
    sil_tlims(test_idx));

%% Plotting
% Plot the example testing reconstruction for one trial
t = (0:length(recon_full{1})-1)/fs; % time array, for plotting the envelopes
figure
hold on
plot(t,stims{test_idx(1)},'k'); % original envelope
plot(t,recon_full{1},'b'); % reconstructed envelope, including silences
recon_nosil_fulllength = NaN(length(recon_full{1}),1); % create an array of NaNs
recon_nosil_fulllength(~sil_chk{test_idx(1)}) = recon_nosil{1}; % insert the reconstruction
    % this is so that silences won't be plotted
plot(t,recon_nosil_fulllength,'r');
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Envelope');
legend('Original','Recon with full envelope','Recon without silences');
% Display the reconstruction accuracies in the title
tle = sprintf('Trail %d: r_{full}=%.3f, r_{no silences}=%.3f',...
    all_trials(test_idx(1)),r_full(1),r_nosil(1));
title(tle);