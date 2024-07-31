function mTRF_plotForwardTRF_UI(chanlocs,models,rAllElec,rAllTrials,figureHandle,plotTypes)
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

    if ~exist('plotTypes') || isempty(plotTypes)
        plotTypes.rAll = 'ScatterPlot';
        plotTypes.normTRF = 1;
    end
    
    set(0, 'currentfigure', figureHandle); 
    rAll = mean(rAllElec,1); % averaging across electrodes
    tmin = models(1).t(1);
    tmax = models(1).t(end);
    
    rAllTrialsVec = rAllTrials{1};
    for iiSub = 2:length(rAllTrials)
        if size(rAllTrialsVec,1) ~= size(rAllTrials{1},1)
            rAllTrialsVec = [];
            break
        end
        rAllTrialsVec(:,iiSub) = rAllTrials{iiSub};
    end
    
    % Plot average TRF
    plotNormFlag = plotTypes.normTRF;
    avgModel = mTRFmodelAvg(models,plotNormFlag);
    
    % Plot EEG prediction correlation
    subplot(2,3,1)
    if strcmp(plotTypes.rAll,'ScatterPlot')
        plot(ones(length(rAll),1),rAll,'.k','MarkerSize',15)
        if min(0,min(rAll)) < max(rAll)
            ylim([min(0,min(rAll)),max(rAll)])
        end
    elseif strcmp(plotTypes.rAll,'BarPlot')
        bar(rAll,'k')
        hold on
        errorbar(mean(rAll),std(rAll),'.','Color','k','LineWidth',3)
    elseif strcmp(plotTypes.rAll,'BoxPlot')
        boxplot(rAll)
        yLim = get(gca,'YLim');
        set(gca,'YLim', [min(yLim(1),0) yLim(2)]);
    end
    xlim([0.5,1.5])
    xticks([])
    ylabel('Neural prediction (r)')
    prepExport()
    grid on
    
    % Plot EEG prediction correlation per trial
    subplot(2,3,2)
    if ~isempty(rAllTrialsVec)
        if strcmp(plotTypes.rAll,'ScatterPlot')
            plot(mean(rAllTrialsVec,2),'.:k','MarkerSize',10)
%             if min(0,min(rAll)) < max(rAll)
%                 ylim([min(0,min(rAll)),max(rAll)])
%             end
        elseif strcmp(plotTypes.rAll,'BarPlot')
%             bar(rAll,'k')
%             hold on
%             errorbar(mean(rAll),std(rAll),'.','Color','k','LineWidth',3)
        elseif strcmp(plotTypes.rAll,'BoxPlot')
%             boxplot(rAll)
%             yLim = get(gca,'YLim');
%             set(gca,'YLim', [min(yLim(1),0) yLim(2)]);
        end
%         xlim([0.5,1.5])
%         xticks([])
%         ylabel('Neural prediction (r)')
%         prepExport()
%         grid on
    end
    
    if ~isempty(chanlocs)
        subplot(2,3,3)
        topoplot(mean(rAllElec,2),chanlocs,'electrodes','off');
        title('Prediction Corr (Avg)')
        caxis([-0.1,0.1])
        colorbar
        prepExport()
    end
    
    % Plot avg TRF model
    subplot(2,3,4)
    plot(avgModel.t,squeeze(avgModel.w))
    title('Env TRF (avg)')
    xlabel('Time-latency (ms)')
    ylabel('Magnitude (a.u.)')
    xlim([tmin+50,tmax-50])
    if plotNormFlag
        ylim([-4,4])
    end
    prepExport()
    grid on

    % Plot GFP
    subplot(2,3,5)
    mTRFplot(avgModel,'gfp',[],'all');
    title('TRF Global Field Power')
    prepExport()
    grid on
end
