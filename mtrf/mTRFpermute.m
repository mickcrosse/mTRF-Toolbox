function null_stats = mTRFpermute(y,pred,method,varargin)
% STATS = MTRFPERMUTE(Y,PRED,METHOD)
% Shuffle or circularly the true signals (y) relative to the predictions
% (pred) in order to calculate a null distribution of accuracies. Note that
% the input format is similar to mTRFevaluate, which calculates the
% correlation or error between the true signal (y) and the prediction
% (pred). Both y and pred should be cell arrays where each cell is a trial
% or fold.
% METHOD has two options:
% - 'permute': Randomly permute the trials or folds
% - 'circshift': Keep trial pairings, but randomly circularly shift the
%       stimuli in each trial or fold
% (Update 30-5-2025: NZ) The inputs are different that what can be found in
% `Crosse et al (2021) Front. Neurosci.` After running some simulations, I
% found that it was not necessary to re-create the model for each
% permutation. Simply permuting the true and predicted signals is enough to
% get the correct null distribution, and it is much more efficient.
% Nate Zuk (2025)

nperm = 1000; % number of permutations

if ~isempty(varargin)
    for n = 2:2:length(varargin)
        if strcmp(varargin{n-1},'nperm')
            nperm = varargin{n};
        end
    end
end

ndim = size(y{1},2); % number of output channels, assumed to be the same for all trials/folds

% Check if the method is correctly specified
if ~strcmp(method,'circshift') && ~strcmp(method,'permute')
    error('Method must be circshift or permute');
end

% Calculate the null distribution
fprintf('Computing the null distribution of accuracies (%d iterations)',nperm);
nullcmp_timer = tic; % keep track of how long it takes to run
r = NaN(nperm,ndim);
err = NaN(nperm,ndim);
for n = 1:nperm
    % display a . every 50 trials
    if mod(n,50)==0, fprintf('.'); end

    if strcmp(method,'permute')
        % randomly select a pair of trials
        sidx = randi(length(y));
        other_idx = setxor(1:length(y),sidx); % remove the index for sidx
        eidx = other_idx(randi(length(pred)-1)); % select another index that is not sidx
    else
        sidx = randi(length(y)); % randomly select a trial...
        eidx = sidx; % ...but use the original trial pairings
    end

    if strcmp(method,'circshift')
        k = randi(length(y{sidx})); % randomly select some amount of time shift, up to the length of y
    else
        k = 0; % no shift
    end

    % Get the true and predicted signal pair
    e = pred{eidx};
    s = y{sidx};
    % match the lenghths of the two signals (important if they are from
    % mismatched trials)
    len = min([size(e,1) size(s,1)]);
    e = e(1:len,:); s = s(1:len,:);
    % apply the circular shift to the true signal
    if k~=0
        s = circshift(s,k);
    end

    % Evaluate the fit between the true and predicted signals
    [r(n,:),err(n,:)] = mTRFevaluate(s,e);
end

null_stats.r = r;
null_stats.err = err;

% insert a new line in the command window when this is completed
fprintf('\n-- Completed @ %.3f s\n',toc(nullcmp_timer));