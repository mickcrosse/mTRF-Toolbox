function plotChannelLocation(chanlocs)
% Handy function plotting the position of the channels on an empty topoplot
% Dependencies: EEGLAB
%
% Author: Giovanni Di Liberto
% Last update: 9 July 2021
%
    figure;topoplot(zeros(length(chanlocs),1),chanlocs,'electrodes','numbers')
end