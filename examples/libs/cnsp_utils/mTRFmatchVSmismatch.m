function stats = mTRFmatchVSmismatch(stim,neural,dir,tmin,tmax,lambdas,windowSize)
%MTRFMATCHVSMISMATCH evaluates the TRF ability to identify correct portions
% of the data with cross-validation. For example, backward models would
% reconstruct the stimulus for unseen data, and that reconstruction would
% then be correlated with the actual stimulus and with a random portion of
% the stimulus (mismatch). Higher correlations on the actual stimulus
% correspond to a correct classification. Cross-correlation is performed
% across trials. The classification is performed on the unseen trial after
% partitioning it in chunks of size WindowSize (e.g., 3 seconds).
%
% Note: the mismatch trial is obtained by using a circularshift of the same
% trial. This is valid under the assumption that there is no repetition
% in that trial. In general, this is more appropriate than other approached
% such as random shuffling of all datapoints, as the circular shift
% preserves the serial correlation in the data.
%
% Acknowledgement:
% The match-vs-mismatch classification metric has been often discussed by
% Prof. Alain de Cheveigne as a particularly useful and flexible evaluation
% method. This code is simply implementing his idea, which he also
% discussed at the CogHear workshop 2022, led by Mounya Elhilali,
% Malcolm Slaney, and Shihab Shamma.
%
% A good reference for reflection on this issue is:
% Geoffrey Brookshire, Nature Human Behaviour, 2022
% https://www.nature.com/articles/s41562-022-01364-0
%
%   MTRF_PLOTFORWARDTRF(MODELS)
%
%       'stim'    -- CND structure containing the stimulus
%       'neural'  -- CND structure containing the neural data
%       'dir'     -- 1: Forward TRF; -1: Backward TRF
%       'tmin'    -- minimum time-lag of the regression window in ms
%       'tmax'    -- maximum time-lag of the regression window in ms
%       'lambdas' -- array with the values to explore for the
%                    regularisation parameter (lambda)
%       'windowSize' -- Size of the window for the match-vs-mismatch task
%                       (dafault: the full duration of each trial)
%
%   Author: Giovanni Di Liberto
%   Last update: 24 June 2022
%   Copyright 2022 Di Liberto Lab, Trinity College Dublin

if ~exist('windowSize') || isempty(windowSize)
    windowSize = 0; % the entire trial window is used
end

% Creating mismatch stimulus
stimMismatch = stim;
for tr = 1:length(stim.data)
    shiftSamples = size(stim.data{tr},1)/2; % circular shifts of half of the length of the trial +- 5% samples
    shiftSamples = round((shiftSamples + 2*(rand-0.5)*0.05*shiftSamples));
    stimMismatch.data{tr} = stimMismatch.data{tr}(circshift(1:size(stimMismatch.data{tr},1),shiftSamples),:);
end

if (0)
    % Just to check that the code works, you can try to shuffle also the stim
    % As a result, the match-vs-mismatch becomes a random comparison
    % i.e., ~50% accuracy
    for tr = 1:length(stim.data)
        shiftSamples = size(stim.data{tr},1)/3; % circular shifts of half of the length of the trial +- 5% samples
        shiftSamples = round((shiftSamples + 2*(rand-0.5)*0.05*shiftSamples));
        stim.data{tr} = stim.data{tr}(circshift(1:size(stim.data{tr},1),shiftSamples),:);
    end
end

% Match-vs-mismatch crossval
stats = mTRFattncrossval(stim.data,stimMismatch.data,neural.data,neural.fs,dir,tmin,tmax,lambdas,'verbose',0,'window',windowSize);

