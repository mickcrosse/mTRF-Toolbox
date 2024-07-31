function [shifted_stim,shft] = randshift_discrete(stim,shiftcols)
% This function is intended for stimuli that are nonzero at a particular 
% set of time ponits (like word or note onsets) and zero otherwise.
% Circularly shift the stimulus values on each trial along the columns
% specified. The size of the shift is in the number of nonzero stimulus values
% and is randomized across trials. On each trial the shift is the same for all columns.
% Note: If the columns have different numbers of nonzero values, the same
% shift size in nonzero stimulus values will produce different shift sizes
% in time.
% Nate Zuk (2022)

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

shft = NaN(length(stim),1); % store how much each stimulus is shifted
shifted_stim = cell(length(stim),1);
for ii = 1:length(stim)
    % for now, copy the stimulus to shifted_stim. The onset values will be
    % overwritten later
    shifted_stim{ii} = stim{ii};
    shft(ii) = randi([1,max_shift]);
    % now circularly shift the discrete values
    for jj = 1:length(shiftcols)
        nz_idx = stim{ii}(:,shiftcols(jj))~=0; % find non-zero indexes
        vals = stim{ii}(nz_idx,shiftcols(jj));
        shftvals = circshift(vals,shft(ii));
        % replace the original non-zero values with the shifted ones
        shifted_stim{ii}(nz_idx,shiftcols(jj)) = shftvals;
    end
end