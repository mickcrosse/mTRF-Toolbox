% Simulate the effect of different types of permutations by simulated a bunch
% of trials with an embedded signal created by a convolution with added
% noise.  Use leave-one-out testing, which is closer to what we normally do.
% After testing, get the null distribution by either:
% 1) shuffling the testing data and using the same model to test
% 2) circularly shifting the testing data, but keeping stimulus-response
% pairs ths same and using the original model
% 3) shuffling both training and testing data to recompute a model (with the
% same amount of regularization) and test
% 4) circularly shifting both the training and testing data to recompute a
% model and test
% Nate Zuk (2020, updated 28-5-2025)
% ** (Updated 2025) Simulate a continuous recording session that is
% segmented into trials.

addpath(genpath('../mtrf')); % the path for the toolbox functions
addpath('../data'); % used to get the example EEG data

Fs = 128; % sampling frequency of the signals
dur = 60*30; % duration of the overall simulation (in s)
trial_dur = 60; % duration of each trial (in s)
freq_range = [1 15]; % frequency range of the signals
snr = -30; % signal to noise ratio in the response (in dB)
lambdas = [0 10.^(0:8)]; % set of ridge regularization parameters to using during mTRFcrossval
%lambdas = 100;
nperm = 5000; % number of times to shuffle the data and get null testing values
quantiles_to_plot = [0.05 0.95];

%%% Generate a TRF model (an exponentially decaying sinusoid with a
%%% frequency of 6 Hz
resp_dur = 350; % duration of the response (in ms)
resp_t = (0:ceil(resp_dur/1000*Fs))/Fs;
resp_frq = 6; % frequency of the response
true_trf = sin(2*pi*resp_frq*resp_t).*exp(-resp_t/(resp_dur/1000/2));

%%% Load an example segment of EEG -- we will use this to create EEG-shaped
%%% noise
exmp_data = load('speech_data');
exmp_eeg = exmp_data.resp(:,85); % use only channel Fz of the example data
EXMP_EEG = fft(exmp_eeg); % get the fourier transform of the EEG
dur_idx = size(exmp_eeg,1); % get the duration of the example data, in indexes
clear exmp_data % delete the example data structure, we no longer need it

%% Simulate the data
%%% Simulate low-frequency data, which requires regularization
%%% There should be an embedded signal that is produced by convolution with
%%% some TRF-like filter, and added noise with the same frequency
%%% distribution.
% stimulus is random bandpass noise
disp('Creating the stimuli and responses...');
s = bandlimited_noise(freq_range,dur,Fs); 
lags = round(resp_t*Fs);
X = lagGen(s,lags);
trfout = X*true_trf';

% generate the noise by randomizing the phases of the EEG segment multiple
% times, and concatenating noise segments to match the overall stimulus
% duration
ns = [];
while size(ns,1)<size(s,1)
    rand_ph = exp(1j*rand(dur_idx,1)*2*pi);
    NS = EXMP_EEG.*rand_ph;
    ns_seg = zscore(real(ifft(NS)));
    % ramp the noise segment to avoid discontinuities betwee segments
    ns_seg = rampstim(ns_seg,Fs,1/freq_range(2));
    ns = [ns; ns_seg];
end
% truncate ns to match the duration of the stimulus
ns = ns(1:size(s,1),:);
eeg = trfout + ns/10^(snr/20)*std(trfout);

% segment the stimuli based on the number of trials
ntr = floor(dur/trial_dur);
seg_idx = round(linspace(1,dur*Fs+1,ntr+1));
stim = cell(ntr,1);
resp = cell(ntr,1);
for n = 1:ntr
    idx = seg_idx(n):seg_idx(n+1)-1;
    stim{n} = s(idx);
    resp{n} = eeg(idx);
end

%% Train, with cross-validation to get the optimum model left-out data
disp('Training on true data...')
true_tm = tic;
pred = cell(ntr,1);
mdl = cell(ntr,1);
opt_idx = NaN(ntr,1);
test_r = NaN(ntr,1);
for n = 1:ntr
    test_tr = n;
    train_trs = setxor(1:ntr,n);
    true_cv_stats = mTRFcrossval(stim(train_trs),resp(train_trs),Fs,1,0,resp_dur,lambdas,'verbose',0);
    % identify the optimal lambda based on the average correlation coefficient
    opt_idx(n) = find(mean(true_cv_stats.r)==max(mean(true_cv_stats.r)),1,'first');
%     fprintf('-- Optimal lambda = %d\n',lambdas(opt_idx));
    mdl{n} = mTRFtrain(stim(train_trs),resp(train_trs),Fs,1,0,resp_dur,lambdas(opt_idx(n)),'verbose',0);
    % save the prediction on each iteration (this is equivalent to
    % re-running the prediction on the same model
    [pred{n},test_stats] = mTRFpredict(stim(test_tr),resp(test_tr),mdl{n},'verbose',0);
    test_r(n) = test_stats.r;
end
fprintf('* Completed training and testing on true data in @ %.3f s\n',toc(true_tm));

%% Plot both models to see how similar they are
figure
hold on
plot(resp_t*1000,true_trf/rms(true_trf),'k','LineWidth',2);
for n = 1:ntr
    plot(mdl{n}.t,mdl{n}.w/rms(mdl{n}.w),'b');
end
set(gca,'FontSize',14);
xlabel('Delay (ms)');
ylabel('RMS normalized model weights');

%% Permutations
%%% Permute 1: Shuffle the testing data to get a null distribution of test
%%% values
disp('** Permute 1 (shuffle testing data)**');
perm_i = tic;
shuff_test_r = NaN(nperm,1);
for n = 1:nperm
    % randomly pick two trials to pair for testing
    % note: this includes correct pairings of trials (which produce high
    % prediction accuracies), so you can get large positive error bars
    % here
    shuff_test_tr = randi(ntr);
    shuff_mdl_tr = randi(ntr); 
    shuff_test_r(n) = corr(pred{shuff_mdl_tr},resp{test_tr});
end
fprintf('* Completed @ %.3f s\n',toc(perm_i));

%%% Permute 2: Randomly select a pairing of trials (this time the correct
%%% pairing) but circularly shift the response
disp('** Permute 2 (circularly shift testing data)**');
perm_ii = tic;
circ_test_r = NaN(nperm,1);
for n = 1:nperm
    circ_pair_tr = randi(ntr);
    rnd_shft = randi(length(resp{circ_pair_tr}));
    shft_resp = circshift(resp{circ_pair_tr},rnd_shft);
    circ_test_r(n) = corr(pred{circ_pair_tr},shft_resp);
end
fprintf('* Completed @ %.3f s\n',toc(perm_ii));

%%% Permute 3: Shuffle both the training and testing data, recompute the
%%% model, and test on the left out (shuffled) data
disp('** Permute 3 (shuffle then re-train model)**');
perm_iii = tic;
use_lambda = lambdas(mode(opt_idx));
shuff_all_r = NaN(nperm,1);
for n = 1:nperm
    shuff_trials = randperm(ntr);
    % randomly select a pair for testing
    test_pair = randi(ntr);
    train_pairs = setxor(1:ntr,test_pair);
    shuff_mdl = mTRFtrain(stim(train_pairs),resp(shuff_trials(train_pairs)),...
        Fs,1,0,resp_dur,use_lambda,'verbose',0);
    [~,shuff_stats] = mTRFpredict(stim(test_pair),resp(shuff_trials(test_pair)),...
        shuff_mdl,'verbose',0);
    shuff_all_r(n) = shuff_stats.r;
end
fprintf('* Completed @ %.3f s\n',toc(perm_iii));

%%% Permute 4: Randomly circularly shift the training and testing data (use
%%% the same shift for all trials), recompute the model, and test on the
%%% left out pair of data
disp('** Permute 4 (circularly shift then re-train model)**');
perm_iv = tic;
circ_all_r = NaN(nperm,1);
for n = 1:nperm
    rnd_shft = randi(dur_idx); % use the length of all trials
    % circularly shift all of the responses
    shft_resp = cell(ntr,1);
    for m = 1:ntr
        shft_resp{m} = circshift(resp{m},rnd_shft);
    end
    % randomly select a pair for testing
    test_pair = randi(ntr);
    train_pairs = setxor(1:ntr,test_pair);
    circ_mdl = mTRFtrain(stim(train_pairs),shft_resp(train_pairs),...
        Fs,1,0,resp_dur,use_lambda,'verbose',0);
    [~,circ_stats] = mTRFpredict(stim(test_pair),shft_resp(test_pair),...
        circ_mdl,'verbose',0);
    circ_all_r(n) = circ_stats.r;
end
fprintf('* Completed @ %.3f s\n',toc(perm_iv));

% Show all of the results
fprintf('\n-- Median [%dth %dth]\n',quantiles_to_plot(1)*100,quantiles_to_plot(2)*100);
fprintf('True r: %.3f [%.3f %.3f]\n',median(test_r),quantile(test_r,...
    quantiles_to_plot(1)),quantile(test_r,quantiles_to_plot(2)));
dpr_shuff_test = (mean(test_r)-mean(shuff_test_r))/sqrt(0.5*(var(test_r)+var(shuff_test_r)));
fprintf('Permute 1 (shuffle testing): %.3f [%.3f %.3f], d-prime = %.3f\n',...
    median(shuff_test_r),quantile(shuff_test_r,quantiles_to_plot(1)),...
    quantile(shuff_test_r,quantiles_to_plot(2)),dpr_shuff_test);
dpr_circ_test = (mean(test_r)-mean(circ_test_r))/sqrt(0.5*(var(test_r)+var(circ_test_r)));
fprintf('Permute 2 (circshift testing): %.3f [%.3f %.3f], d-prime = %.3f\n',...
    median(circ_test_r),quantile(circ_test_r,quantiles_to_plot(1)),...
    quantile(circ_test_r,quantiles_to_plot(2)),dpr_circ_test);
dpr_shuff_all = (mean(test_r)-mean(shuff_all_r))/sqrt(0.5*(var(test_r)+var(shuff_all_r)));
fprintf('Permute 3 (shuffle all): %.3f [%.3f %.3f], d-prime = %.3f\n',...
    median(shuff_all_r),quantile(shuff_all_r,quantiles_to_plot(1)),...
    quantile(shuff_all_r,quantiles_to_plot(2)),dpr_shuff_all);
dpr_circ_all = (mean(test_r)-mean(circ_all_r))/sqrt(0.5*(var(test_r)+var(circ_all_r)));
fprintf('Permute 4 (circshift all): %.3f [%.3f %.3f], d-prime = %.3f\n',...
    median(circ_all_r),quantile(circ_all_r,quantiles_to_plot(1)),...
    quantile(circ_all_r,quantiles_to_plot(2)),dpr_circ_all);

% Plot the results
figure
hold on
plot([1 5],[0 0],'k--');
% True
md = median(test_r);
uq = quantile(test_r,quantiles_to_plot(2));
lq = quantile(test_r,quantiles_to_plot(1));
errorbar(1,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 1
md = median(shuff_test_r);
uq = quantile(shuff_test_r,quantiles_to_plot(2));
lq = quantile(shuff_test_r,quantiles_to_plot(1));
errorbar(2,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 2
md = median(circ_test_r);
uq = quantile(circ_test_r,quantiles_to_plot(2));
lq = quantile(circ_test_r,quantiles_to_plot(1));
errorbar(3,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 3
md = median(shuff_all_r);
uq = quantile(shuff_all_r,quantiles_to_plot(2));
lq = quantile(shuff_all_r,quantiles_to_plot(1));
errorbar(4,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);
% Permute 4
md = median(circ_all_r);
uq = quantile(circ_all_r,quantiles_to_plot(2));
lq = quantile(circ_all_r,quantiles_to_plot(1));
errorbar(5,md,md-lq,uq-md,'k.','MarkerSize',20,'LineWidth',2);

set(gca,'FontSize',14,'XTick',1:5,'XTickLabelRotation',45,...
    'XTickLabel',{'true','shuffle test','circshift test','shuffle all','circshift all'});
ylabel('Pearsons r');

%%% Additional functions %%%
function y = bandlimited_noise(modfreq,dur,Fs)
% Taken from https://github.com/natezuk/EEG-analysis-tools.
% Creates a band-limited (brick filter) random signal, so there is no
% energy outside of the frequency range specified
rndenv = randn(dur*Fs,1);
RNDENV = fft(rndenv);
f = (0:length(RNDENV)-1)/length(RNDENV)*Fs;
idx = (f>=modfreq(1)&f<=modfreq(2))|(f<=Fs-modfreq(1)&f>=Fs-modfreq(2));
    % need to filter evenly on both sides of Fs/2 for real signals
RNDENV(~idx) = 0; % remove those frequencies
y = real(ifft(RNDENV));
end