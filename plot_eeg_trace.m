function plot_eeg_trace(EEG,time_or_fs)
% PLOT_EEG_TRACE(EEG,TIME)
% PLOT_EEG_TRACE(EEG,FS)
% Plot the EEG trace, with each channel separated along the y axis
% Inputs:
% - EEG = the EEG (time x channels)
% There are two options for the second input:
% - time = a time array that must be equal in size to the number of time 
%     indexes in EEG
% - fs = sampling rate of the EEG (Hz)
% Nate Zuk (2019)

% Check if using the time array or fs
if length(time_or_fs)>1 % if time_or_fs is an array of numbers
    if size(EEG,1)~=length(time_or_fs) % if the time array doesn't equal the 
       % number of time samples in EEG...
       error('The time array must have the same number of time samples as the EEG');
    else
       t = time_or_fs;
    end
else % if fs is specified, make the time array, assume starting at 0
    t = (0:size(EEG,2)-1)/fs;
end

% Z-score all of the data in EEG (so spacing is more consistent
% irrespective of the variance in the EEG)
mEEG = mean(reshape(EEG,[numel(EEG),1])); % mean
sEEG = std(reshape(EEG,[numel(EEG),1])); % standard deviation
zEEG = (EEG-mEEG)/sEEG; % z-scored eeg
% scale the magnitude of the EEG traces, so that the channels are more
% separable in the plot
eeg_scaling = 1/2;
zEEG = zEEG*eeg_scaling; % divide the magnitude by 

% Plot the EEG trace
figure
hold on
for n = 1:size(EEG,2)
    y_offset = n; % amount of offset along the y axis relative to 0
    plot(t,zEEG(:,n)+y_offset,'k');
end
set(gca,'FontSize',14,'YLim',[-1 size(EEG,2)+1]);
if size(EEG,2)>10 % if there are more than 10 channels...
    % ...label every 5 on the y axis
    set(gca,'YTick',5:5:size(EEG,2));
else
    set(gca,'YTick',1:size(EEG,2));
end
xlabel('Time');
ylabel('EEG channel');