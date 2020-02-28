function mTRFdemos(demo)
%MTRFDEMOS  mTRF-Toolbox demos.
%   MTRFDEMOS(DEMO) runs one of several short demo scripts using an example
%   dataset to demonstrate the usage of the different functions in
%   mTRF-Toolbox. DEMO is a scalar integer specifying which demo to run
%   and includes the following:
%       '1'     VESPA estimation (checkerboard contrast)
%       '2'     TRF estimation (speech envelope)
%       '3'     STRF estimation (speech spectrogram)
%       '4'     Stimulus reconstruction (speech envelope)
%       '5'     Backward to forward model transformation (speech envelope)
%       '6'     Sinigle-lag stimulus reconstruction (speech envelope)
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

switch demo
    
    
    case 1 % VESPA estimation (checkerboard contrast)
        
        
        % Load data
        load('data/contrast_data.mat','stim','resp','fs','Nf','factor');
        
        % Scale data achieve model weights in uV
        stim = stim*Nf; % checkerboard frame rate (Nf): 60Hz
        resp = resp*factor; % BioSemi ADC factor: 524.288mV / 2^24bits
        
        % Model parameters
        dir = 1;
        tmin = -150;
        tmax = 450;
        lambda = 0; % no regularization needed for stochastic signal
        
        % Compute model weights
        model = mTRFtrain(stim,resp,fs,dir,tmin,tmax,lambda,'zeropad',0);
        
        % Define ROI
        chan = 23; % channel Oz
        
        % Plot TRF & GFP
        figure(1)
        subplot(1,2,1)
        plot(model.t,squeeze(model.w(:,:,chan)),'linewidth',3)
        xlim([-50,350])
        title('Contrast TRF (Oz)')
        xlabel('Time lag (ms)')
        ylabel('Amplitude (\muV)')
        axis square
        grid on
        subplot(1,2,2)
        area(model.t,squeeze(std(model.w,[],3)),'edgecolor','none');
        xlim([-50,350])
        title('Global Field Power')
        xlabel('Time lag (ms)')
        axis square
        grid on
        
        
    case 2 % TRF estimation (speech envelope)
        
        
        % Load data
        load('data/speech_data.mat','stim','resp','fs','factor');
        
        % Normalize data
        stim = sum(stim,2);
        resp = resp*factor;
        
        % Model parameters
        dir = 1;
        tmin = -150;
        tmax = 450;
        lambda = 0.1;
        
        % Compute model weights
        model = mTRFtrain(stim,resp,fs,dir,tmin,tmax,lambda,...
            'method','Tikhonov','zeropad',0);
        
        % Define ROI
        chan = 85; % channel Fz
        
        % Plot TRF & GFP
        figure(2)
        subplot(1,2,1)
        plot(model.t,squeeze(model.w(:,:,chan)),'linewidth',3)
        xlim([-50,350]);
        title('Speech TRF (Fz)')
        xlabel('Time lag (ms)')
        ylabel('Amplitude (a.u.)')
        axis square
        grid on
        subplot(1,2,2)
        area(model.t,squeeze(std(model.w,[],3)),'edgecolor','none');
        xlim([-50,350])
        title('Global Field Power')
        xlabel('Time lag (ms)')
        axis square
        grid on
        
        
    case 3 % STRF estimation (speech spectrogram)
        
        
        % Load data
        load('data/speech_data.mat','stim','resp','fs','factor');
        
        % Normalize data
        resp = resp*factor;
        n = size(stim,2);
        
        % Model parameters
        dir = 1;
        tmin = -150;
        tmax = 450;
        lambda = 0.5;
        
        % Compute model weights
        model = mTRFtrain(stim,resp,fs,dir,tmin,tmax,lambda,...
            'method','ridge','zeropad',0);
        
        % Define ROI
        chan = 85; % channel Fz
        
        % Plot STRF & GFP
        figure(3)
        subplot(2,2,1)
        imagesc(model.t(14:66),1:16,squeeze(model.w(:,14:66,chan)))
        xlim([-50,350])
        set(gca,'ydir','normal')
        title('Speech STRF (Fz)')
        ylabel('Frequency band')
        axis square
        subplot(2,2,2)
        imagesc(model.t(14:66),1:16,squeeze(std(model.w(:,14:66,:),[],3)))
        xlim([-50,350])
        set(gca,'ydir','normal')
        title('Global Field Power')
        axis square
        
        % Compute broadband TRF by summing over STRF channels:
        
        % TRF parameters
        dir = 1;
        tmin = -150;
        tmax = 450;
        lambda = 0.05;
        
        % Compute model weights
        model = mTRFtrain(stim,resp,fs,dir,tmin,tmax,lambda,...
            'method','ridge','zeropad',0);
        
        % Define ROI
        chan = 85; % channel Fz
        
        % Plot TRF & GFP
        subplot(2,2,3)
        plot(model.t,squeeze(sum(model.w(:,:,chan))),'linewidth',3)
        xlim([-50,350])
        title('Speech TRF (Fz)')
        xlabel('Time lag (ms)')
        ylabel('Amplitude (a.u.)')
        axis square
        grid on
        subplot(2,2,4)
        area(model.t,squeeze(std(sum(model.w),[],3)),'edgecolor','none')
        xlim([-50,350])
        title('Global Field Power')
        xlabel('Time lag (ms)')
        axis square
        grid on
        
        
    case 4 % Stimulus reconstruction (speech envelope)
        
        
        % Load data
        load('data/speech_data.mat','stim','resp','fs');
        
        % Normalize data
        stim = sum(stim,2);
        stim = stim/std(stim);
        resp = resp/std(resp(:));
        
        % Downsample data
        fsNew = 64;
        stim = resample(stim,fsNew,fs);
        resp = resample(resp,fsNew,fs);
        fs = fsNew;
        
        % ---Create training/test sets---
        
        %TODO(mickcrosse): create cvfold() function
        
        % Determine folds for cross-validation
        nfold = 10;
        batch = ceil(size(stim,1)/nfold);
        
        % Training set
        stimtrain = cell(nfold,1);
        resptrain = cell(nfold,1);
        for i = 1:nfold
            idx = batch*(i-1)+1:min(batch*i,size(resp,1));
            stimtrain{i} = stim(idx,:);
            resptrain{i} = resp(idx,:);
        end
        
        % Test set
        itest = 1;
        stimtest = stimtrain{itest};
        resptest = resptrain{itest};
        stimtrain(itest) = []; % remove test set from training data
        resptrain(itest) = [];
        
        % ---Cross-validation---
        
        % Cross-validation parameters
        dir = -1;
        tmin = 0;
        tmax = 250;
        lambdas = 10.^(-6:2:6);
        nlambda = length(lambdas);
        
        % Run cross-validation
        cv = mTRFcrossval(stimtrain,resptrain,fs,dir,tmin,tmax,lambdas,...
            'zeropad',0,'fast',1);
        
        % ---Model training---
        
        % Model parameters
        [rmax,idx] = max(mean(cv.acc));
        lambda = lambdas(idx);
        
        % Train model
        model = mTRFtrain(stimtrain,resptrain,fs,dir,tmin,tmax,lambda,...
            'zeropad',0);
        
        % ---Model testing---
        
        % Test model
        [pred,test] = mTRFpredict(stimtest,resptest,model,'zeropad',0);
        
        % Plot cross-validation results
        figure(4)
        subplot(2,2,1)
        errorbar(1:nlambda,mean(cv.acc),std(cv.acc)/sqrt(nfold-1),...
            'linewidth',2)
        set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6)
        xlim([0,nlambda+1])
        title('CV Accuracy')
        xlabel('Lambda (1\times10^\lambda)')
        ylabel('Correlation')
        axis square
        grid on
        subplot(2,2,2)
        errorbar(1:nlambda,mean(cv.err),std(cv.err)/sqrt(nfold-1),...
            'linewidth',2)
        set(gca,'xtick',1:nlambda,'xticklabel',-6:2:6)
        xlim([0,nlambda+1])
        title('CV Error')
        xlabel('Lambda (1\times10^\lambda)')
        ylabel('MSE')
        axis square
        grid on
        
        % Plot test results
        subplot(2,2,3)
        plot(stimtest,'linewidth',2)
        hold on
        plot(pred,'linewidth',2)
        hold off
        xlim([0,400])
        title('Reconstruction')
        xlabel('Sample number')
        ylabel('Amplitude (a.u.)')
        axis square
        grid on
        legend('Orig','Pred')
        subplot(2,2,4)
        bar(1,rmax)
        hold on
        bar(2,test.acc)
        title('Test Result')
        xlabel('Metric')
        ylabel('Correlation')
        set(gca,'xtick',1:2,'xticklabel',{'CV','Test'})
        hold off
        axis square
        grid on
        
        
    case 5 % Backward to forward model transformation (speech envelope)
        
        
        % Load data
        load('data/speech_data.mat','stim','resp','fs');
        
        % Normalize data
        stim = sum(stim,2);
        stim = stim/std(stim);
        resp = resp/std(resp(:));
        
        % Model parameters
        dir = -1;
        tmin = -150;
        tmax = 450;
        lambda = 50;
        
        % Compute backward model weights
        bmodel = mTRFtrain(stim,resp,fs,dir,tmin,tmax,lambda);
        
        % Transform to forward model weights
        fmodel = mTRFtransform(bmodel,resp);
        
        % Define ROI
        chan = 85; % channel Fz
        
        % Plot TRF & GFP
        figure(5)
        subplot(1,2,1)
        plot(fmodel.t,squeeze(fmodel.w(chan,:)),'linewidth',3);
        xlim([-50,350]);
        title('Frontal TRF (Fz)')
        xlabel('Time lag (ms)')
        ylabel('Amplitude (a.u.)')
        axis square
        grid on
        subplot(1,2,2)
        area(fmodel.t,squeeze(std(fmodel.w)),'edgecolor','none');
        xlim([-50,350])
        title('Global Field Power')
        xlabel('Time lag (ms)')
        axis square
        grid on
        
        
    case 6 % Sinigle-lag stimulus reconstruction (speech envelope)
        
        
        % Load data
        load('data/speech_data.mat','stim','resp','fs');
        
        % Normalize data
        stim = sum(stim,2);
        stim = stim/std(stim(:));
        resp = resp/std(resp(:));
        
        % Downsample data
        fsNew = 64;
        stim = resample(stim,fsNew,fs);
        resp = resample(resp,fsNew,fs);
        fs = fsNew;
        
        % ---Create training/test sets---
        
        % Determine folds for cross-validation
        nfold = 10;
        batch = ceil(size(stim,1)/nfold);
        
        % Training set
        stimtrain = cell(nfold,1);
        resptrain = cell(nfold,1);
        for i = 1:nfold
            idx = batch*(i-1)+1:min(batch*i,size(resp,1));
            stimtrain{i} = stim(idx,:);
            resptrain{i} = resp(idx,:);
        end
        
        % ---Cross-validation---
        
        % Cross-validation parameters
        dir = -1;
        tmin = 0;
        tmax = 1000;
        lambdas = 10.^-2;
        
        % Run cross-validation
        [stats,t] = mTRFcrossval(stimtrain,resptrain,fs,dir,tmin,tmax,...
            lambdas,'type','single','zeropad',0);
        
        mr = squeeze(mean(stats.acc))';
        vr = squeeze(var(stats.acc))';
        mrmse = squeeze(mean(stats.err))';
        vrmse = squeeze(var(stats.err))';
        
        % Plot reconstruction accuracy
        figure(6)
        subplot(1,2,1)
        error = [fliplr(mr-sqrt(vr/nfold)),mr+sqrt(vr/nfold)];
        x = [-fliplr(t),-t];
        h = fill(x,error,'b','edgecolor','none');
        set(h,'facealpha',0.2)
        hold on
        plot(-fliplr(t),fliplr(mr),'linewidth',2)
        hold off
        xlim([tmin,tmax])
        title('Reconstruction Accuracy')
        xlabel('Time lag (ms)')
        ylabel('Correlation')
        axis square
        grid on
        subplot(1,2,2)
        error = [fliplr(mrmse-sqrt(vrmse/nfold)),mrmse+sqrt(vrmse/nfold)];
        x = [-fliplr(t),-t];
        h = fill(x,error,'b','edgecolor','none');
        set(h,'facealpha',0.2)
        hold on
        plot(-fliplr(t),fliplr(mrmse),'linewidth',2)
        hold off
        xlim([tmin,tmax])
        title('Reconstruction Error')
        xlabel('Time lag (ms)')
        ylabel('MSE')
        axis square
        grid on
        
        
end