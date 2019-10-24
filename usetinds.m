function tinds = usetinds(tlims,Fs,maxind)
% tinds = usetinds(tlims,Fs,maxind)
% Identify the indexes that should be used based on the time limits
% provided, the sampling rate, and the stimulus duration. It is assumed
% that the starting index is t=0. The output is a logical array with 1 for
% each index to be used.
% Inputs:
% - tlims = specifies the range of time to use (in s):
%       tlims=[] --> all indexes are used
%       tlims=[t(1)] --> t(1) to <end_of_array>-t(1) are
%           included
%       tlims=[t(1) Inf] --> t(1) to end of array
%       tlims=[0 t(1)] --> beginning of array to t(1)
%       tlims=[t(1) t(2)] --> t(1) to t(2)
%       tlims=t(1,...,n>2) --> samples at the specific time points
%           specified in tlims
% - Fs = sampling rate (Hz)
% - maxind = maximum index to include (should be set to the length of
%       data being sampled)
% Outputs:
% - tinds = logical array with maxind indexes with 1 for each index that
%       falls within the range/set of times specified in tlims
% Nate Zuk (2017)

t = (0:maxind-1)/Fs;
if ~isempty(tlims),
    if length(tlims)==1,
        % set limits relative to start and end of stimulus
        tinds = t>=tlims&t<=(maxind/Fs-tlims);
    elseif length(tlims)==2,
        % use the times in usets as limits
        tinds = t>=tlims(1)&t<=tlims(2);
    elseif length(tlims)>2,
        tinds = false(1,maxind);
        tidx = round(tlims*Fs);
        tidx(tidx>maxind) = []; % remove indexes that are greater than maxind
        tinds(tidx) = true;
    end
else
    tinds = true(1,maxind);
end