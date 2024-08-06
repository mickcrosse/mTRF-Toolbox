function mTRF_plotForwardTRF(neural,models,rAllElec,accMM)
%MTRF_PLOTBACKWARDTRF It plots standard results, such as EEG prediction
% correlations and mTRF weights. The figure must have been created already.
%   MTRF_PLOTBACKWARDTRF(MODELS)
%
%       'neural' -- CND structure containing the neural data (the field
%                   chanlocs will be used
%       'models' -- structure with the mTRF models (output of mTRFtrain)
%                   for multiple participants
%       'rAllElec' -- array containing the average prediction correlations
%                   for each subject and electrode.
%       'accMM'  -- accuracy of a match-vs-mismatch task. Values from 0 to
%                   1 (optional)
%
%   Author: Giovanni Di Liberto
%   Last update: 23 June 2022
%   Copyright 2022 Di Liberto Lab, Trinity College Dublin

    if exist('accMM') && ~isempty(accMM)
        nSubplots = 3;
    else
        nSubplots = 2;
    end
    
    rAll = mean(rAllElec,1); % averaging across electrodes
    tmin = models(1).t(1);
    tmax = models(1).t(end);
    
    % Plot average TRF
    plotNormFlag = 1;
    avgModel = mTRFmodelAvg(models,plotNormFlag);
    
    % Plot EEG prediction correlation
    subplot(1,nSubplots,1)
    plot(ones(length(rAll),1),rAll,'.k','MarkerSize',15)
    xlim([0.5,1.5])
    xticks([])
%     ylim([min(0,min(rAll)),max(rAll)])
    ylabel('Reconstruction corr (r)')
    run prepExport.m
    grid on 
       
    % Plot avg TRF model
    subplot(1,nSubplots,nSubplots)
    plot(avgModel.t,squeeze(avgModel.w))
    title('Model weights')
    xlabel('Time-latency (ms)')
    ylabel('Magnitude (a.u.)')
    xlim([tmin+50,tmax-50])
    ylim([-4,4])
    run prepExport.m
    grid on
    
    if nSubplots > 2
        % Plot avg TRF model
        subplot(1,nSubplots,2)
        plot(ones(length(accMM),1),accMM*100,'.k','MarkerSize',15)
        xlim([0.5,1.5])
        xticks([])
        ylim([0,100])
%         ylim([min(0,min(accMM)),max(accMM)])
        ylabel('MM classification (%)')
        run prepExport.m
        grid on 
    end
end
