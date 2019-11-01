function plot_trf(plot_title,forward_trf,fs,minlag,maxlag,chan_select)
% PLOT_TRF(TRF,FS,MINLAG,MAXLAG,CHAN_SELECT)
% Plot a forward TRF model in the following ways:
%  If there is only one feature:
%  - plots the EEG trace of all selected channels separated along the y-axis 
%    (traditional EEG trace plot)
%  - a "butterfly" plot showing all selected EEG channels on the same y axis
%  - the average across all selected channels
%  - the global field power across all selected channels (standard
%    deviation across channels)
% ** Make sure the model is a forward model (predicting EEG). If you are
%    using a forward model output by mTRFtransform, set constant_term_check
%    to 0
% Inputs:
% - plot_title = title for all of the plots
% - forward_trf = the forward TRF model (lags x channels)
% - fs = sampling rate (Hz)
% - minlag = minimum lag in the TRF model (ms)
% - maxlag = maximum lag in the TRF model (ms)
% - chan_select (optional) = list of channels to plot
%       default: plot all channels
% Nate Zuk (2019)

% if channels aren't specified, use all of them
if nargin<6, chan_select = 1:size(forward_trf,2); end

% identifies if the model contains a constant term
% (it does by default after running mTRFtrain or mTRFcrossval, but
% mTRFtransform does not retain a constant term)
% if nargin<7, constant_term_check = true;  end

% Create the lags array
lags = (floor(minlag/1000*fs):ceil(maxlag/1000*fs))/fs*1000;

% Reshape the trf into channels x lags
% if constant_term_check % if the first index of each channel is the constant term
% %     TRF = reshape(forward_trf,[nchans length(lags)+1]);
% %     C = forward_trf(1,:);
%     TRF = forward_trf(2:end,:);
% else
%     TRF = reshape(trf,[nchans length(lags)]);
TRF = squeeze(forward_trf); % dimension 1 is the input, which will always be one for a
    % univariate TRF (like envelope)
%     C = [];
% end

% Plot the EEG trace
plot_eeg_trace(TRF(:,chan_select),lags);
% set(gca,'YTickLabel',chan_select); % label the y axis with the correct channel numbers
xlabel('Lag (ms)');
title(plot_title)

% Plot the butterfly plot
figure
plot(lags,TRF(:,chan_select));
set(gca,'FontSize',16);
xlabel('Lag (ms)');
ylabel('TRF (A.U.)');
title(plot_title)

% Plot the EEG averaged across channels
figure
subplot(1,2,1);
plot(lags,mean(TRF(:,chan_select),2),'k','LineWidth',2);
set(gca,'FontSize',16);
xlabel('Lag (ms)');
ylabel('Average across channels (A.U.)');
title(plot_title)

% Plot the GFP
subplot(1,2,2);
plot(lags,std(TRF(:,chan_select),[],2),'k','LineWidth',2);
set(gca,'FontSize',16);
xlabel('Lag (ms)');
ylabel('Global field power (A.U.)');
title(plot_title)