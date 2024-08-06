function [nChans,nChansVec] = cndGetNChans(cnd,verbose)
%CNDGETNCHANS  Get the number of channels of the neural signal in a CND
%              structure.
%   NCHANS = CNDGETNCHANS(CND) returns the number of channels of the neural
%   signal in the CND structure. The function also checks that the number
%   of channels is consistent across participants. In that case, the value
%   -1 is assigned to 'nChans' and a warning is displayed indicating to
%   look into the variable 'nChansVec', which indicates the number of
%   channels for each subject in the dataset.
%
%   Author: Giovanni Di Liberto
%   Last update: 11 June 2021
%   Copyright 2021 Di Liberto Lab, Trinity College Dublin

    if ~exist('verbose')
        verbose = 1;
    end
    
    if isempty(cnd) || isempty(cnd.data)
        disp('The CND structure is empty or not a cell array')
        return
    elseif ~iscell(cnd.data)
        disp('The CND.data structure is not a cell array')
        return
    end
    
    nChans = size(cnd.data{1},2);
    
    nSubs = length(cnd.data);
    nChansVec = zeros(nSubs,1);
    for sub = 1:nSubs
        nChansVec(sub) = size(cnd.data{sub},2);
    end
    if verbose && sum(nChansVec ~= nChans)
        disp('Warning: Subjects have a different number of channels. Please check the variable "nChansVec"')
        nChans = -1;
    end
end