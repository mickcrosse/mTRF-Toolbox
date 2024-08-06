function eeg_interped = removeBadChannels(eeg, chanlocs, rejectedChannels,highThreshold,lowThreshold)
% This function identifies and replaces bad channels with a spline
% interpolation of all other channels.
% Dependencies: EEGLAB.
%
% Input:
%  eeg           - eeg data (time x channel)
%  chanlocs      - channel location strucuture
%  rejectedChannels - if empty, channels to reject are calculated based on
%                     the variance of the eeg data; else, only the channels
%                     specified in this vector are rejected
%  highThreshold - this function rejects channels with variance larger than
%                  VAR*highThreshold, where VAR is the median variance
%                  across all channels (default 3)
%                  
%  lowThreshold  - this function rejects channels with variance larger than
%                  VAR/lowThreshold, where VAR is the median variance
%                  across all channels (default 6 - very conservative)
%
% Output:
%  eeg_interped  - eeg after bad channel interpolation
%
% Example:
%  eeg = removeBadChannels(eeg, chanlocs, []);
%  eeg = removeBadChannels(eeg, chanlocs, [], 4);
%  eeg = removeBadChannels(eeg, chanlocs, [], 4, 10);
%
% Originally by Edmund Lalor - 04/09/09.
% Edited by Giovanni Di Liberto - 12/05/14

    if ~exist('highThreshold') || isempty(highThreshold)
        highThreshold = 3;
    end
    
    if ~exist('lowThreshold') || isempty(lowThreshold)
        lowThreshold = 6; % conservative threshold
    end
    
    eeg = eeg';
    
    EEG.data = eeg;
    EEG.chanlocs = chanlocs;
    EEG.nbchan = length(EEG.chanlocs);
    EEG.setname = 'EEGLAB structure';
    EEG.pnts = size(eeg,2);
    EEG.xmax = size(eeg,2);
    EEG.xmin = 1;

    EEG.trials = 1;
    EEG.srate = 1;
    EEG.chaninfo = [];
    EEG.icawinv = [];
    EEG.icaweights = [];
    EEG.icasphere  = [];
    EEG.icaact = [];
    EEG.epoch = [];
    EEG.specdata = [];
    EEG.icachansind = [];
    EEG.specicaact = [];
    EEG.reject = [];
    EEG.stats = [];
    EEG.ref = 'averef';
    EEG.etc = [];
    EEG.event = [];

    std_chans = zeros(1,size(eeg,1));
    for i = 1:size(eeg,1) 
        std_chans(i) = std(eeg(i,:));
    end
    mean(std_chans);

    if ~exist('rejectedChannels','var') || isempty(rejectedChannels)
        % Identify bad channels
        badchans = [];
        
        for i = 1:size(eeg,1)
            % If the std of the channel is more than 3 times the mean of
            % the standard deviations of all channels
            if std(eeg(i,:)) > mean(std_chans)*highThreshold 
                badchans = [badchans i];
            end
        end
        
        clear std_chans
        for i = 1:size(eeg,1)
            if ~isempty(find(badchans==i, 1))
                continue
            end
            std_chans(i) = std(eeg(i,:));
        end
        mean(std_chans);
        
        % If the standard deviation of the channel is smaller than
        % 1/lowThreshold times the mean of the stds of all channels
        for i = 1:size(eeg,1)
            if std(eeg(i,:)) < mean(std_chans)/lowThreshold 
                badchans = [badchans i];
            end
        end
    else
        badchans = rejectedChannels;
    end

    if (~isempty(badchans))
        iEEG = eeg_interp(EEG, badchans); % check that eeglab is included in the dependencies in order to run this function
        disp(['badChannels: ' num2str(badchans)])
        eeg_interped = iEEG.data;
        eeg_interped = eeg_interped(1:size(eeg,1),:);
    else
        eeg_interped = eeg;
    end

    for i = 1:size(eeg,1) 
        std_chans(i) = std(eeg_interped(i,:));
    end

%     run_std = mean(std_chans);
    eeg_interped = eeg_interped';
end

