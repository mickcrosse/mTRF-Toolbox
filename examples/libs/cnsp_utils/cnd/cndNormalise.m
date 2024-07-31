function [cnd] = cndNormalise(cnd,celldim)
%CNDNORMALISE Rescales stimulus or neural data in CND format. The function
% will implement a number of normalisations and standardisations. This
% first version only implements one type of standardisation, consisting of
% dividing the data by the overall standard deviation i.e., std calculated
% across all dimensions (e.g., EEG channels). This procedure preserves the
% relative magnitude across feature dimensions (e.g., freq bands of a
% speech sgram, EEG channels).
%   CND = CNDNORMALISE(CND,CELLDIM) 
%
%      'cnd'     -- stimulus or neural data in the Continuous-event Neural
%                    Data format (CND)
%      'celldim' -- dimensions to normalise. If stimulus data, dimensions
%                   correspond to different stimulus feature-sets
%                   (e.g., sgram, word dissimilarity)
%
%   Note: cnd.data: featureSets x trials. 'celldim' refers to the first
%   dimension of that cell array.
%
%   Author: Giovanni Di Liberto
%   Last update: 23 June 2022
%   Copyright 2022 Di Liberto Lab, Trinity College Dublin

    if isempty(cnd) || isempty(cnd.data)
        disp('The CND structure is empty or not a cell array')
        return
    elseif ~iscell(cnd.data)
        disp('The CND.data structure is not a cell array')
        return
    end
    
    % If celldim was unspecified or empty
    if ~exist('celldim') || isempty(celldim)
        celldim = 1:size(cnd.data,1);
    end
    
    % For each celldim (e.g., stimulus feature)
    for iCell = celldim
        % Calculating normalisation factor
        tmpCnd = cnd.data{iCell,:};
        for tr = 2:length(cnd.data) % getting all values
            tmpCnd = cat(1,tmpCnd,cnd.data{tr});
        end
        tmpCnd = tmpCnd(:);
        normFactor = std(tmpCnd);

        % Warning if the data appears to have too many zeros, as this
        % standardisation may not be what we want to apply
        if sum(tmpCnd==0)/length(tmpCnd) > 0.75
            disp(['Warning: dimension ',num2str(iCell),' has over 3 times more zeros than non-zero values. Skipping normalisation.'])
        else
            % Rescaling
            for tr = 1:length(cnd.data)
                cnd.data{iCell,tr} = cnd.data{iCell,tr}/normFactor;
            end
        end
        clear tmpCnd;
    end
end
