function [pred,stats] = mTRFpredict(stim,resp,model,varargin)
%MTRFPREDICT  Predict and evaluate the output of a model.
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
%       'Dir'       -- direction of causality (forward=1, backward=-1)
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
%   [PRED,STATS] = MTRFPREDICT(...) returns the test statistics in a
%   structure with the following fields:
%       'r'         -- prediction accuracy based on Pearson's correlation
%                      coefficient (nfold-by-yvar)
%       'err'       -- prediction error based on the mean squared error
%                      (nfold-by-yvar)
%
%   [...] = MTRFPREDICT(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both STIM and RESP.
%       'corr'      A string specifying the correlation metric to use:
%                       'Pearson'   Pearson's linear correlation
%                                   coefficient (default): suitable for
%                                   data with a linear relationship
%                       'Spearman'  Spearman's rank correlation
%                                   coefficient: suitable for data with a
%                                   non-linear relationship
%       'error'     A string specifying the error metric to use:
%                       'mse'       mean square error (default): take the
%                                   square root to convert it to the
%                                   original units of the data (i.e., RMSE)
%                       'mae'       mean absolute error: more robust to
%                                   outliers than MSE
%       'split'     A scalar specifying the number of segments in which to
%                   split each trial of data when computing the covariance
%                   matrices. This is useful for reducing memory usage on
%                   large datasets. By default, the entire trial is used.
%       'window'    A scalar specifying the window size over which to
%                   compute performance in seconds. By default, the entire
%                   trial or segment is used.
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%       'verbose'   A numeric or logical specifying whether to execute in
%                   verbose mode: pass in 1 for verbose mode (default), or
%                   0 for non-verbose mode.
%
%   See also PREDICT, CORR, MSE, MAE, MTRFCROSSVAL.
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

%   Authors: Mick Crosse <crossemj@tcd.ie>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Edmund Lalor <edlalor@tcd.ie>
%            Nate Zuk <zukn@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Define X and Y variables
if model.Dir == 1
    x = stim; y = resp;
elseif model.Dir == -1
    x = resp; y = stim;
end

% Format data in cell arrays
[x,xobs,xvar] = formatcells(x,arg.dim,arg.split);
[y,yobs,yvar] = formatcells(y,arg.dim,arg.split);
if isempty(y)
    yobs = xobs;
    yvar = size(model.w,3);
end

% Check equal number of observations
if ~isequal(xobs,yobs)
    error(['STIM and RESP arguments must have the same number of '...
        'observations.'])
end

% Convert time lags to samples
tmin = model.t(1)*model.fs/1e3;
tmax = model.t(end)*model.fs/1e3;
lags = tmin:tmax;
arg.window = round(arg.window*model.fs);

% Compute sampling interval
delta = 1/model.fs;

% Get dimensions
xvar = unique(xvar);
yvar = unique(yvar);
nfold = numel(x);
switch model.type
    case 'multi'
        nlag = 1;
    case 'single'
        nlag = numel(lags);
end

% Truncate output
if ~arg.zeropad && ~isempty(y)
    [y,yobs] = truncate(y,tmin,tmax,yobs);
end
if arg.window
    nwin = sum(floor(yobs/arg.window));
else
    nwin = nfold;
end

% Verbose mode
if arg.verbose
    v = verbosemode([],sum(yobs),0,nfold);
end

% Format model weights
switch model.type
    case 'multi'
        model.w = [model.b;reshape(model.w,[xvar*numel(lags),yvar])]*delta;
    case 'single'
        model.w = [model.b;model.w]*delta;
end

% Initialize variables
pred = cell(nfold,1);
if nargout > 1 && ~isempty(y)
    r = zeros(nwin,yvar,nlag);
    err = zeros(nwin,yvar,nlag);
end
ii = 0;

% Test model
for i = 1:nfold
    
    if arg.window
        ii = ii(end)+1:ii(end)+floor(yobs(i)/arg.window);
    else
        ii = i;
    end
    
    % Test set
    xlag = lagGen(x{i},lags,arg.zeropad,1);
    
    switch model.type
        
        case 'multi'
            
            % Predict output
            pred{i} = xlag*model.w;
            
            if nargout > 1 && ~isempty(y)
                
                % Evaluate performance
                [r(ii,:),err(ii,:)] = mTRFevaluate(y{i},pred{i},...
                    'corr',arg.corr,'error',arg.error,...
                    'window',arg.window);
                
            end
            
        case 'single'
            
            % Initialize cell
            pred{i} = zeros(yobs(i),yvar,nlag);
            
            for j = 1:nlag
                
                % Index lag
                idx = [1,xvar*(j-1)+2:xvar*j+1];
                
                % Predict output
                pred{i}(:,:,j) = xlag(:,idx)*squeeze(model.w(:,j,:));
                
                if nargout > 1 && ~isempty(y)
                    
                    % Evaluate performance
                    [r(ii,:,j),err(ii,:,j)] = ...
                        mTRFevaluate(y{i},pred{i}(:,:,j),...
                        'corr',arg.corr,'error',arg.error,...
                        'window',arg.window);
                    
                end
                
            end
            
    end
    
    % Verbose mode
    if arg.verbose
        v = verbosemode(v,[],i,nfold);
    end
    
end

% Format output arguments
if nfold == 1
    pred = pred{1};
end
if nargout > 1 && ~isempty(y)
    
    stats = struct('r',r,'err',err);
    
    % Verbose mode
    if arg.verbose
        verbosemode(v,[],i+1,nfold,stats);
    end
    
end

function v = verbosemode(v,nobs,fold,nfold,stats)
%VERBOSEMODE  Execute verbose mode.
%   V = VERBOSEMODE(V,NOBS,FOLD,NFOLD,STATS) prints details about the
%   progress of the main function to the screen.

if fold == 0
    v = struct('msg',[],'h',[],'tocs',0);
    fprintf('\nTest on %d samples\n',nobs)
    fprintf('Predicting output\n')
    v.msg = ['%d/%d [',repmat(' ',1,nfold),']\n'];
    v.h = fprintf(v.msg,fold,nfold);
elseif fold <= nfold
    if fold == 1 && toc < 0.1
        pause(0.1)
    end
    v.tocs = v.tocs + toc;
    fprintf(repmat('\b',1,v.h))
    v.msg = ['%d/%d [',repmat('=',1,fold),repmat(' ',1,nfold-fold),'] - ',...
        '%.3fs/fold\n'];
    v.h = fprintf(v.msg,fold,nfold,v.tocs/fold);
end
if fold < nfold
    tic
elseif fold > nfold
    rmax = mean(stats.r,1); rmax = max(rmax(:));
    emax = mean(stats.err,1); emax = max(emax(:));
    fprintf('correlation: %.4f - error: %.4f\n',rmax,emax)
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

% Correlation metric
corrOptions = {'Pearson','Spearman'};
validFcn = @(x) any(validatestring(x,corrOptions));
addParameter(p,'corr','Pearson',validFcn);

% Error metric
errOptions = {'mse','mae'};
validFcn = @(x) any(validatestring(x,errOptions));
addParameter(p,'error','mse',validFcn);

% Split data
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'split',1,validFcn);

% Window size
errorMsg = 'It must be a positive numeric scalar within indexing range.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'window',0,validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'zeropad',true,validFcn); % zero-pad design matrix
addParameter(p,'verbose',true,validFcn); % verbose mode

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially-matched strings
arg.corr = validatestring(arg.corr,corrOptions);
arg.error = validatestring(arg.error,errOptions);