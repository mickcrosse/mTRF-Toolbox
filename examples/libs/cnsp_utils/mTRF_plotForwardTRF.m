function mTRF_plotForwardTRF(neural,models,rAllElec)
%MTRF_PLOTFORWARDTRF It plots standard results, such as EEG prediction
% correlations and mTRF weights. The figure must have been created already.
%   MTRF_PLOTFORWARDTRF(MODELS)
%
%       'neural' -- CND structure containing the neural data (the field
%                   chanlocs will be used
%       'models' -- structure with the mTRF models (output of mTRFtrain)
%                   for multiple participants
%       'rAllElec'   -- array containing the average prediction correlations
%                   for each subject and electrode.
%   Author: Giovanni Di Liberto
%   Last update: 23 June 2022
%   Copyright 2022 Di Liberto Lab, Trinity College Dublin

    rAll = mean(rAllElec,1); % averaging across electrodes
    tmin = models(1).t(1);
    tmax = models(1).t(end);
    
    % Plot average TRF
    plotNormFlag = 1;
    avgModel = mTRFmodelAvg(models,plotNormFlag);
    
    % Plot EEG prediction correlation
    subplot(2,2,1)
    plot(ones(length(rAll),1),rAll,'.k','MarkerSize',15)
    xlim([0.5,1.5])
    xticks([])
    if min(0,min(rAll)) < max(rAll)
        ylim([min(0,min(rAll)),max(rAll)])
    end
    ylabel('Prediction corr (r)')
    run prepExport.m
    grid on
    
    if isfield(neural,'chanlocs')
        subplot(2,2,2)
        title('Prediction Corr (Avg)')
        topoplot(mean(rAllElec,2),neural.chanlocs,'electrodes','off');
        caxis([-0.2,0.2])
        colorbar
        run prepExport.m
    end
    
    % Plot avg TRF model
    subplot(2,2,3)
    plot(avgModel.t,squeeze(avgModel.w))
    title('Env TRF (avg)')
    xlabel('Time-latency (ms)')
    ylabel('Magnitude (a.u.)')
    xlim([tmin+50,tmax-50])
    ylim([-4,4])
    run prepExport.m
    grid on

    % Plot GFP
    subplot(2,2,4)
    mTRFplot(avgModel,'gfp',[],'all');
    title('TRF Global Field Power')
    run prepExport.m
    grid on
end
