function stats = mTRFpermute(stim,resp,fs,Dir,method,tmin,tmax,lambda,varargin)
% STATS = MTRFPERMUTE(STIM,RESP,FS,DIR,METHOD,TMIN,TMAX,LAMBDA)
% Shuffle or circularly shift stimuli relative to the response and
% recalculate a null distribution of r values for a given lambda value
% When running this, include all of the folds or trials in 'stim' and
% 'resp'.
% METHOD has three options:
% - 'permute': Randomly permute the trials or folds
% - 'circshift': Keep trial pairings, but randomly circularly shift the
%       stimuli in each trial or fold
% - 'both': Both randomly permute trials or folds and randomly circularly 
%       shift the stimuli
% Nate Zuk (2022)

nperm = 100; % number of permutations

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

% Check if the method is correctly specified
if ~strcmp(method,'circshift') && ~strcmp(method,'permute') && ~strcmp(method,'both')
    error('Method must be circshift, permute, or both');
end

fprintf('Computing the null distribution of accuracies (%d iterations)\n',nperm);
nullcmp_timer = tic; % keep track of how long it takes to run
nullr = NaN(nperm,ndim);
nullerr = NaN(nperm,ndim);
for n = 1:nperm
    % display a . every 10 trials
    if mod(n,10)==0, fprintf('(%d/%d)\n',n,nperm); end

    if strcmp(method,'permute') || strcmp(method,'both')
        % randomly select pairs of trials
        eidx = randperm(length(resp));
        sidx = randperm(length(stim));
    else
        eidx = randperm(length(resp)); % randomly shuffle...
        sidx = eidx; % ...but use the original trial pairings
    end
    if strcmp(method,'circshift') || strcmp(method,'both')
        k = NaN(length(stim),1);
        for ii = 1:length(stim)
            % randomly select a circular shift amount
            shift_range = [ceil(tmax/1000*fs) size(stim{sidx(ii)},1)+floor(tmin/1000*fs)];
            k(ii) = randi(shift_range);
        end
    else
        k = zeros(length(stim),1); % no shift
    end

    % Create the stim and resp arrays for this iteration
    e = resp(eidx);
    s = stim(sidx);
    for ii = 1:length(stim)
        len = min([size(e{ii},1) size(s{ii},1)]);
        e{ii} = e{ii}(1:len,:); s{ii} = s{ii}(1:len,:);
        if k(ii)~=0
            s{ii} = circshift(s{ii},k(ii));
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

end

stats.nullr = nullr;
stats.nullerr = nullerr;

% insert a new line in the command window when this is completed
fprintf('-- Completed @ %.3f s\n',toc(nullcmp_timer));