% Simulate data of stimulus tracking (1-15 Hz) with added EEG-shaped noise
% for 20 2-minute-long trials, calculate a null distribution of prediction
% accuracies, and compare the true predictiona accuracies to the null
% prediction accuracies.
% Nate Zuk (2025)

addpath(genpath('../mtrf')); % the path for the toolbox functions
addpath('../data'); % used to get the example EEG data

%% Simulation parameters
Fs = 128; % sampling frequency of the signals
ntr = 20; % number of trials
freq_range = [1 15]; % frequency range of the signals
snr = -30; % signal to noise ratio in the response (in dB)
lambdas = [0 10.^(0:8)]; % set of ridge regularization parameters to using during mTRFcrossval
nperm = 500; % number of times to shuffle the data and get null testing values

%% Generate the TRF model
%%% Generate a TRF model (an exponentially decaying sinusoid with a
%%% frequency of 6 Hz
resp_dur = 350; % duration of the response (in ms)
resp_t = (0:ceil(resp_dur/1000*Fs))/Fs;
resp_frq = 6; % frequency of the response
true_trf = sin(2*pi*resp_frq*resp_t).*exp(-resp_t/(resp_dur/1000/2));

%% Load example EEG for EEG-shaped noise
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
stim = cell(ntr,1);
resp = cell(ntr,1);
lags = round(resp_t*Fs);
for n = 1:ntr
    % create a stimulus that has equal energy within the frequency range
    % specified (by default this is 1-15 Hz)
    s = bandlimited_noise(freq_range,dur_idx/Fs,Fs); 
    % set the variance of the original stimulus to 1
    stim{n} = zscore(s);
    % Convolve the stimulus with the true TRF
    X = lagGen(stim{n},lags);
    trfout = X*true_trf';
    % add EEG-shaped noise
    rand_ph = exp(1j*rand(dur_idx,1)*2*pi);
    NS = EXMP_EEG.*rand_ph; % phase shift the signal based on the random values
    ns = zscore(real(ifft(NS))); % convert back into time domain and normalize
    resp{n} = trfout+ns/10^(snr/20)*std(trfout); % scale the noise to the appropriate SNR
end

%% Iteratively leave one trial out, cross-validate on the rest
%%% Train, with cross-validation to get the optimum model, and test on
%%% left-out data
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
disp('Creating the null distribution');
null_stats = mTRFpermute(resp,pred,'circshift');

%% Plot the true values and the null distribution
figure
% plot r values
hold on
[h_r,edges] = histcounts(null_stats.r,20); % create the histogram with 20 bins
r_bins = edges(1:end-1)+diff(edges)/2;
bar(r_bins,h_r/sum(h_r),1,'k');
stem(test_r,0.4*ones(length(test_r),1),'Color','r','LineWidth',2);
xlabel('Prediction accuracy (r)');
ylabel('Prop. of null distribution');
legend('Null distribution','True values','Location','northwest');


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