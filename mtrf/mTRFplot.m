function h = mTRFplot(model,type,feat,chan,xlims,avgfeat,avgchan)
%MTRFPLOT  Plot TRF model weights.
%   H = MTRFPLOT(MODEL) plots the weights of a TRF model as a function of
%   time lag. By default, this functon sums over all features and plots
%   all channels individually. To modify the proerties of the plot, use
%   set(h,'PARAM1',VAL1,'PARAM2',VAL2,...).
%
%   H = MTRFPLOT(MODEL,TYPE) specifies the type of plot to generate. Pass
%   in 'trf' for TYPE to generate a standard TRF plot, or 'gfp' to generate
%   an area plot of the global field power (GFP). Pass in 'mtrf' or 'mgfp'
%   to generate an image plot of a multivariate TRF or GFP, repectively.
%
%   H = MTRFPLOT(MODEL,TYPE,FEAT,CHAN) specifies the feature and channel
%   numbers to plot, respectively. Pass in a vector of integer values for
%   FEAT or CHAN to specify them, or 'all' to use all of them. By default,
%   all features and channels are used.
%
%   H = MTRFPLOT(MODEL,TYPE,FEAT,CHAN,XLIMS) specifies the x-axis limits.
%   Setting this within this function results in better scaling of the
%   multivariate image plots.
%
%   H = MTRFPLOT(MODEL,TYPE,FEAT,CHAN,XLIMS,AVGFEAT,AVGCHAN) specifies
%   whether to average over the selected features and channels,
%   respectively. Pass in 1 for AVGFEAT or AVGCHAN to average them, or 0 to
%   plot them individually.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Set default values
if nargin < 2 || isempty(type)
    type = 'trf';
end
if model.Dir == -1
    model.w = permute(model.w,[3,2,1]);
end
if nargin < 3 || isempty(feat) || strcmpi(feat,'all')
    feat = 1:size(model.w,1);
end
if nargin < 4 || isempty(chan) || strcmpi(chan,'all')
    chan = 1:size(model.w,3);
end
if nargin < 5 || isempty(xlims)
    xlims = [model.t(1),model.t(end)];
end
if nargin < 6 || isempty(avgfeat)
    avgfeat = 1;
end
if nargin < 7 || isempty(avgchan)
    avgchan = 0;
end

% Define lags
switch type
    case {'mtrf','mgfp'}
        [~,idx1] = min(abs(model.t-xlims(1)));
        [~,idx2] = min(abs(model.t-xlims(2)));
        lags = idx1:idx2;
end

% Define features and channels
model.w = model.w(feat,:,chan);

% Average features
switch type
    case {'trf','gfp'}
        if avgfeat
            model.w = mean(model.w,1);
        end
end

% Average channels
switch type
    case {'trf','mtrf'}
        if avgchan
            model.w = mean(model.w,3);
        end
    case {'gfp','mgfp'}
        model.w = std(model.w,1,3);
end

% Generate plot
switch type
    case 'trf'
        h = plot(model.t,squeeze(model.w),'linewidth',3); grid on
    case 'gfp'
        h = area(model.t,model.w,'edgecolor','none'); grid on
    case {'mtrf','mgfp'}
        h = imagesc(model.t(lags),1:numel(feat),model.w(:,lags,:));
        set(gca,'ydir','normal')
end
xlabel('Time lag (ms)')
xlim(xlims)
axis square