function [TRF,C] = plot_multifeature_trf(plot_title,forward_trf,fs,minlag,maxlag,feats,chan_select,constant_term_check)
% PLOT_TRF(TRF,FS,MINLAG,MAXLAG,CHAN_SELECT)
% Plot a forward TRF model for multiple features:
%  - In two color images, plot the average across EEG channels and the
%    global field power (standard deviation across channels), with each
%    feature separated along the y axis
%  - Plot the average and GFP in separate plots, using the same y axis
% ** Make sure the model is a forward model (predicting EEG). If you are
%    using a forward model output by mTRFtransform, set constant_term_check
%    to 0
% Inputs:
% - plot_title = title for all of the plots
% - trf = the forward TRF model (lags x channels)
% - fs = sampling rate (Hz)
% - minlag = minimum lag in the TRF model (ms)
% - maxlag = maximum lag in the TRF model (ms)
% - feats = a list of feature names or values (can be a cell array of strings 
%       or an integer array)
% - chan_select = list of channels to plot
%       default: plot all channels
% - constant_term_check = 1 if the model contains a constant term, 0
%     otherwise (default = 1; set this to 0 if the model was output by
%     mTRFtransform)
% Outputs:
% - TRF = reformatted TRF (channels x lags)
% - C = array of constant terms (one for each feature)
% Nate Zuk (2019)

% if channels aren't specified, use all of them
if nargin<7, chan_select = 1:size(forward_trf,2); end

% identifies if the model contains a constant term
% (it does by default after running mTRFtrain or mTRFcrossval, but
% mTRFtransform does not retain a constant term)
if nargin<8, constant_term_check = true;  end

% Create the lags array
lags = (floor(minlag/1000*fs):ceil(maxlag/1000*fs))/fs*1000;

% Reshape the trf into channels x lags
if constant_term_check % if the first index of each channel is the constant term
    C = forward_trf(1,:);
    TRF = reshape(forward_trf(2:end,:),[length(feats) length(lags) size(forward_trf,2)]);
else % only use if the forward model came from mTRFtransform
    TRF = reshape(forward_trf,[length(feats) length(lags) size(forward_trf,2)]);
    C = [];
end

% Plot an image of the average across channels, where each channel is
% separated along the y axis
figure
subplot(1,2,1);
imagesc(lags,1:length(feats),squeeze(mean(TRF(:,:,chan_select),3)));
axis('xy')
colorbar;
set(gca,'FontSize',16,'YTick',1:length(feats),'YTickLabel',feats);
xlabel('Lag (ms)');
ylabel('Features');
title([plot_title ', average across channels']);

% Do the same for the GFP
subplot(1,2,2);
imagesc(lags,1:length(feats),squeeze(std(TRF(:,:,chan_select),[],3)));
axis('xy')
colorbar;
set(gca,'FontSize',16,'YTick',1:length(feats),'YTickLabel',feats);
xlabel('Lag (ms)');
ylabel('Features');
title([plot_title ', GFP']);

% Plot the TRF averaged across channels, with each feature overlaid on the same axis
% make the legend for these plots
if isnumeric(feats)
    leg_feats = cell(length(feats),1);
    for n = 1:length(feats)
        leg_feats{n} = num2str(feats(n));
    end
else
    leg_feats = feats;
end
figure
cmap = colormap('jet');
subplot(1,2,1);
hold on
for n = 1:length(feats)
    clr_idx = round((n-1)/length(feats)*size(cmap,1))+1; % index of color to use in colormap
        % scaled so that colors are evenly distributed across all possible
        % colors in the map
    plot(lags,squeeze(mean(TRF(n,:,chan_select),3)),'Color',cmap(clr_idx,:),...
        'LineWidth',2);
end
set(gca,'FontSize',16);
xlabel('Lag (ms)');
ylabel('Average across channels (A.U.)');
legend(leg_feats);
title(plot_title)

% Plot the GFP
subplot(1,2,2);
hold on
for n = 1:length(feats)
    clr_idx = round((n-1)/length(feats)*size(cmap,1))+1; % index of color to use in colormap
        % scaled so that colors are evenly distributed across all possible
        % colors in the map
    plot(lags,squeeze(std(TRF(n,:,chan_select),[],3)),'Color',cmap(clr_idx,:),...
        'LineWidth',2);
end
set(gca,'FontSize',16);
xlabel('Lag (ms)');
ylabel('Average across channels (A.U.)');
legend(leg_feats);
title(plot_title)