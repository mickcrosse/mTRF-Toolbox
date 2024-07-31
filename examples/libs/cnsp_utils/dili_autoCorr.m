function dili_autoCorr(x,fs)
% This function plots the autocorrelation of a signal x with dimension
% (timeSamples x 1).
%
% Author: Giovanni Di Liberto
% Last update: 9 July 2021
%
    xACorr = xcorr(x,x);
    xACorr = xACorr(1:round(length(xACorr+1)/2));
    xACorr = xACorr(end-fs*20:end);
    figure;plot(flip(-(0:length(xACorr)-1)/fs),xACorr)
end