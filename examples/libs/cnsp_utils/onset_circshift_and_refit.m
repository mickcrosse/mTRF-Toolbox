function stats = onset_circshift_and_refit(stim,resp,fs,Dir,shiftcols,tmin,tmax,lambda,varargin)
% This function is intended for stimuli that are nonzero at a particular 
% set of time ponits (like word or note onsets) and zero otherwise.
% Circularly shift the stimulus values (in # nonzero indexes) and
% recalculate a null distribution of r values for a given lambda value
% When running this, include all of the folds or trials in 'stim' and
% 'resp'.
% Nate Zuk (2022)

nperm = 100; % number of permutations

% Parse varargin (can use used to modify the nperm parameter
if ~isempty(varargin)
    for n = 2:2:length(varargin)
        if strcmp(varargin{n-1},'nperm')
            nperm = varargin{n};
        end
    end
end

if Dir==1
    ndim = size(resp{1},2); % number of output channels, assumed to be the same for all trials
elseif Dir==-1
    ndim = size(stim{1},2);
end

% Identify the maximum size of the circular shift, based on the minimum
% number of non-zero values of the columns that should be shifted
ntr = length(stim); % number of trials
nnonzero = NaN(length(shiftcols),ntr);
for ii = 1:ntr
    for jj = 1:length(shiftcols)
        nnonzero(jj,ii) = sum(stim{ii}(:,shiftcols(jj))~=0);
    end
end
max_shift = min(min(nnonzero));

fprintf('Computing the null distribution of accuracies (%d iterations)\n',nperm);
nullcmp_timer = tic; % keep track of how long it takes to run
nullr = NaN(nperm,ndim);
nullerr = NaN(nperm,ndim);
testtr = NaN(nperm,1);
for n = 1:nperm
    % display progress every 10 trials
    if mod(n,10)==0, fprintf('(%d/%d)\n',n,nperm); end

    % Randomly reorder trials -- the last one will be used for testing
    eidx = randperm(length(resp)); % randomly shuffle...
    sidx = eidx; % ...but use the original trial pairings

    % Identify the amount of shift for each trial
    k = NaN(length(stim),1);
    for ii = 1:length(stim)
        % randomly select a circular shift amount
        shift_range = [1 max_shift];
        k(ii) = randi(shift_range);
    end

    % Create the stim and resp arrays for this iteration
    e = resp(eidx);
    s = stim(sidx);
    for ii = 1:length(stim)
        %circularly shift the discrete values
        for jj = 1:length(shiftcols)
            nz_idx = s{ii}(:,shiftcols(jj))~=0; % find non-zero indexes
            vals = s{ii}(nz_idx,shiftcols(jj));
            shftvals = circshift(vals,k(ii));
            % replace the original non-zero values with the shifted ones
            s{ii}(nz_idx,shiftcols(jj)) = shftvals;
        end
    end

    % Compute the model on all trials with one left out
    tridx = 1:length(s)-1;
    tstidx = length(s);
    mdl = mTRFtrain(s(tridx),e(tridx),fs,Dir,tmin,tmax,lambda,'verbose',0);
    % Test the model on the left out trial
    [~,iter_stats] = mTRFpredict(s{tstidx},e{tstidx},mdl,'verbose',0);
    nullr(n,:) = iter_stats.r;
    nullerr(n,:) = iter_stats.err;
    testtr(n) = sidx(end); % get the testing trial used for this iteration

end

% Save the results
stats.nullr = nullr;
stats.nullerr = nullerr;
stats.testtr = testtr;

% insert a new line in the command window when this is completed
fprintf('-- Completed @ %.3f s\n',toc(nullcmp_timer));