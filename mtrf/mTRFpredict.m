function [pred,stats] = mTRFpredict(stim,resp,model,varargin)
%MTRFPREDICT  Predict model output.
%   PRED = MTRFPREDICT(STIM,RESP,MODEL) predicts the output of a forward
%   encoding model (stimulus to neural response) or a backward decoding
%   model (neural response to stimulus) using time-lagged input features.
%   STIM and RESP are matrices or cell arrays containing corresponding
%   trials of continuous test data.
%
%   MODEL is a structure with the following fields:
%       'w'         -- normalized model weights (xvar-by-nlag-by-yvar)
%       'b'         -- normalized bias term (1-by-nlag-by-yvar)
%       't'         -- time lags (ms)
%       'fs'        -- sample rate (Hz)
%       'dir'       -- direction of causality (forward=1, backward=-1)
%       'type'      -- type of model (multi-lag, single-lag)
%
%   If STIM or RESP are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   STIM and RESP must have the same number of observations.
%
%   If STIM or RESP are cell arrays containing multiple trials, separate
%   predictions are made for each trial and returned in a cell array. STIM
%   and RESP must contain the same number of trials.
%
%   [PRED,STATS] = MTRFPREDICT(...) returns the statistics in a structure
%   with the following fields:
%       'acc'       -- prediction accuracy based on Pearson's correlation
%                      coefficient (ntrial-by-yvar)
%       'err'       -- prediction error based on the mean squared error
%                      (ntrial-by-yvar)
%
%   [...] = MTRFPREDICT(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both STIM and RESP.
%       'acc'       A string specifying the accuracy metric to use:
%                       'Pearson'   Pearson's linear correlation
%                                   coefficient (default): suitable for
%                                   data with a linear relationship
%                       'Spearman'  Spearman's rank correlation
%                                   coefficient: suitable for data with
%                                   non-linear relationship
%       'err'       A string specifying the error metric to use:
%                       'msc'       Mean square error (default): take the
%                                   square root to convert it to the
%                                   original units of the data (i.e., RMSE)
%                       'mae'       Mean absolute error: more robust to
%                                   outliers than MSE
%       'split'     A scalar specifying the number of segments in which to
%                   split each trial of data when computing the covariance
%                   matrices. This is useful for reducing memory usage on
%                   large datasets. By default, the entire trial is used.
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%
%   See mTRFdemos for examples of use.
%
%   See also MTRFTRAIN, MTRFTRANSFORM, MTRFCROSSVAL, MTRFMULTICROSSVAL.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Lalor EC (2009) Modeling the human visual system using the
%          white-noise approach. 2009 4th International IEEE/EMBS
%          Conference on Neural Engineering, Antalya 589-592.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Edmund Lalor <edmundlalor@gmail.com>
%            Nate Zuk <zukn@tcd.ie>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Define X and Y variables
if model.dir == 1
    x = stim; y = resp;
elseif model.dir == -1
    x = resp; y = stim;
end

% Format data in cell arrays
[x,xobs,xvar] = formatcells(x,arg.dim);
[y,yobs,yvar] = formatcells(y,arg.dim);
if isempty(y)
    yobs = xobs;
    yvar = size(model.w,3);
end

% Check equal number of observations
if ~isequal(xobs,yobs)
    error('STIM and RESP arguments must have the same number of observations.')
end

% Convert time lags to samples
tmin = model.t(1)*model.fs/1e3;
tmax = model.t(end)*model.fs/1e3;
lags = tmin:tmax;

% Compute sampling interval
delta = 1/model.fs;

% Get dimensions
nlag = numel(lags);
xvar = unique(xvar);
yvar = unique(yvar);
ntrial = numel(x);

% Format model weights
switch model.type
    case 'multi'
        model.w = [model.b;reshape(model.w,[xvar*nlag,yvar])]*delta;
        nlag = 1;
    case 'single'
        model.w = [model.b;model.w]*delta;
end

% Initialize performance variables
pred = cell(ntrial*arg.split,1);
acc = zeros(ntrial*arg.split,yvar,nlag);
err = zeros(ntrial*arg.split,yvar,nlag);

% Test model
n = 0;
for i = 1:ntrial
    
    % Max segment size
    nseg = ceil(xobs(i)/arg.split);
    
    for j = 1:arg.split
        
        n = n+1;
        
        % Test data
        iseg = nseg*(j-1)+1:min(nseg*j,xobs(i));
        [xlag,idx] = lagGen(x{i}(iseg,:),lags,arg.zeropad);
        xlag = [ones(numel(idx),1),xlag]; %#ok<*AGROW>
        
        switch model.type
            
            case 'multi'
                
                % Predict output
                pred{n} = xlag*model.w;
                
                % Evaluate performance
                if nargout > 1
                    [acc(n,:),err(n,:)] = mTRFevaluate(y{i}(idx,:),...
                        pred{n},'acc',arg.acc,'err',arg.err);
                end
                
            case 'single'
                
                pred{n} = zeros(numel(idx),yvar,nlag);
                
                for k = 1:nlag
                    
                    % Predict output
                    ilag = [1,xvar*(k-1)+2:xvar*k+1];
                    pred{n}(:,:,k) = xlag(:,ilag)*squeeze(model.w(:,k,:));
                    
                    % Evaluate performance
                    if nargout > 1
                        [acc(n,:,k),err(n,:,k)] = ...
                            mTRFevaluate(y{i}(idx,:),pred{n}(:,:,k),...
                            'acc',arg.acc,'err',arg.err);
                    end
                    
                end
                
        end
        
    end
    
end

% Format output arguments
if ntrial == 1
    pred = pred{:};
end
if nargout > 1
    stats = struct('acc',acc,'err',err);
end

function arg = parsevarargin(varargin)
%PARSEVARARGIN  Parse input arguments.
%   [PARAM1,PARAM2,...] = PARSEVARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   parses the input arguments of the main function.

% Create parser object
p = inputParser;

% Dimension to work along
errorMsg = 'It must be a positive integer scalar within indexing range.';
validFcn = @(x) assert(x==1||x==2,errorMsg);
addParameter(p,'dim',1,validFcn);

% Accuracy metric
accOptions = {'Pearson','Spearman'};
validFcn = @(x) any(validatestring(x,accOptions));
addParameter(p,'acc','Pearson',validFcn);

% Error metric
errOptions = {'mse','mae'};
validFcn = @(x) any(validatestring(x,errOptions));
addParameter(p,'err','mse',validFcn);

% Split data
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'split',1,validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'zeropad',true,validFcn); % zero-pad design matrix
addParameter(p,'gpu',false,validFcn); % run on GPU

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially-matched strings
arg.acc = validatestring(arg.acc,accOptions);
arg.err = validatestring(arg.err,errOptions);