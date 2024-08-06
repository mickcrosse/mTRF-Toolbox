function mTRFplotTuningCurve(stats,lambdas,dir)
% This function plots the tuning curve of the mTRF cross-validation
% procedure.
%
% Input:
%  stats   - output of mTRFcrossval containing the evaluation metric
%  lambdas - the lambda values used as an input of mTRFcrossval
%  dir     - direction of the mTRF model; input of mTRFcrossval
%
% Author: Giovanni Di Liberto
% Last update: 9 July 2021
%
    figure;
    if dir>0
        plot(mean(mean(stats.r,1),3),'.-k','MarkerSize',10);
    else
        plot(mean(stats.r),'.-k','MarkerSize',10);
    end
    xlabel('Lambda')
    xticks(1:2:length(lambdas))
    xticklabels(lambdas(1:2:end))
    if dir>0
        ylabel('Neural response prediction correlation (r)')
    else
        ylabel('Stimulus reconstruction correlation (r)')
    end
    title('Tuning curve')
    run prepExport.m
    drawnow;
end