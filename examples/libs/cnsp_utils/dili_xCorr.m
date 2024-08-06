function [xXCorr,timeVec] = dili_xCorr(x,y,fs,winSec)
% This function plots the cross-correlation between two signals x and y
% (both with size timeSamples x 1).
%
% Author: Giovanni Di Liberto
% Last update: 9 July 2021
%
    xXCorr = xcorr(x,y,'coeff');
    middle = round(length(xXCorr+1)/2);
%     xXCorr = xXCorr(1:round(length(xXCorr+1)/2));
    xXCorr = xXCorr(middle-fs*winSec:middle+fs*winSec);
    timeVec = (-ceil(length(xXCorr)/2)+1:ceil(length(xXCorr)/2)-1)/fs;
    figure;plot(timeVec,xXCorr)
end