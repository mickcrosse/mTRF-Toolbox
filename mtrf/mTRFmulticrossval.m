function [stats,t] = mTRFmulticrossval(stim,resp,resp1,resp2,fs,Dir,tmin,tmax,lambda,varargin)
%MTRFMULTICROSSVAL  Cross-validation for multisensory model optimization.
%   STATS = MTRFMULTICROSSVAL(STIM,RESP,RESP1,RESP2,FS,DIR,TMIN,TMAX,LAMBDA)
%   cross validates a forward encoding model (stimulus to neural response)
%   or a backward decoding model (neural response to stimulus) over
%   multiple trials of data for optimizing an additive model of
%   multisensory processing. Additive models are trained on the summed
%   covariances of the unisensory responses RESP1 and RESP2, and validated
%   on the actual multisensory responses RESP as per Crosse et al. (2015).
%   Pass in 1 for DIR to validate a forward model, or -1 to validate a
%   backward model. STIM, RESP, RESP1 and RESP2 are cell arrays containing
%   corresponding trials of continuous data. FS is a scalar specifying the
%   sample rate in Hertz, and TMIN  and TMAX are scalars specifying the
%   minimum and maximum time lags in milliseconds. For backward models, the
%   time lags are automatically reversed. LAMBDA is a vector of
%   regularization values to be validated and controls overfitting.
%
%   MTRFMULTICROSSVAL returns the cross-validation statistics in a
%   structure with the following fields:
%       'r'         -- correlation coefficient based on Pearson's linear
%                      correlation coefficient (nfold-by-nlambda-by-yvar)
%       'err'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar)
%
%   MTRFMULTICROSSVAL performs a leave-one-out cross-validation over all
%   trials. To achieve a k-fold cross-validation, arrange STIM, RESP, RESP1
%   and RESP2 in k-by-1 cell arrays. The number of folds can also be
%   increased by an integer factor using the 'split' parameter (see below).
%
%   If STIM, RESP, RESP1 or RESP2 contain matrices, it is assumed that the
%   rows correspond to observations and the columns to variables, unless
%   otherwise stated via the 'dim' parameter (see below). If they contain
%   vectors, it is assumed that the first non-singleton dimension
%   corresponds to observations. Each trial of STIM, RESP, RESP1 or RESP2
%   must have the same number of observations.
%
%   [STATS,T] = MTRFMULTICROSSVAL(...) returns a vector containing the time
%   lags used in milliseconds. These data are useful for interpreting the
%   results of single-lag models.
%
%   [...] = MTRFMULTICROSSVAL(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
%   specifies additional parameters and their values. Valid parameters are
%   the following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both STIM and RESP.
%       'method'    A string specifying the regularization method to use:
%                       'ridge'     ridge regression (default): suitable
%                                   for multivariate input features
%                       'Tikhonov'  Tikhonov regularization: dampens fast
%                                   oscillatory components of the solution
%                                   but may cause cross-channel leakage for
%                                   multivariate input features
%                       'ols'       ordinary least squares: equivalent to
%                                   setting LAMBDA=0 (no regularization)
%       'type'      A string specifying the type of model to fit:
%                       'multi'     use all lags simultaneously to fit a
%                                   multi-lag model (default)
%                       'single'    use each lag individually to fit
%                                   separate single-lag models
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
%       'fast'      A numeric or logical specifying whether to use the fast
%                   cross-validation method (requires more memory) or the
%                   slower method (requires less memory): pass in 1 to use
%                   the fast method (default), or 0 to use the slower
%                   method. Note, both methods are numerically equivalent.
%       'verbose'   A numeric or logical specifying whether to execute in
%                   verbose mode: pass in 1 for verbose mode (default), or
%                   0 for non-verbose mode.
%
%   Notes:
%   Each iteration of MTRFMULTICROSSVAL partitions the N trials or segments
%   of data into two subsets, fitting a model to N-1 trials (training set)
%   and validating it on the left-out trial (validation set). Performance
%   on the validation set can be used to optimize hyperparameters (e.g.,
%   LAMBDA, TMAX). Once completed, it is recommended to evaluate model
%   performance on separate held-out data using the mTRFpredict function.
%
%   Discontinuous trials of data should not be concatenated prior to cross-
%   validation, as this will introduce artifacts in places where the
%   temporal integration window crosses over trial boundaries. Each trial
%   of continuous data should be input as a separate cell.
%
%   See also CROSSVAL, MTRFPARTITION, MTRFMULTITRAIN, MTRFPREDICT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Butler JS, Lalor EC (2015) Congruent Visual Speech
%          Enhances Cortical Entrainment to Continuous Auditory Speech in
%          Noise-Free Conditions. J Neurosci 35(42):14195-14204.
%      [2] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <crossemj@tcd.ie>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Edmund Lalor <edlalor@tcd.ie>
%            Nate Zuk <zukn@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Validate input parameters
validateparamin(fs,Dir,tmin,tmax,lambda)

% Define X and Y variables
if Dir == 1
    x = stim; y = resp;
elseif Dir == -1
    x = resp; y = stim;
    [tmin,tmax] = deal(tmax,tmin);
end

% Format data in cell arrays
[x,xobs,xvar] = formatcells(x,arg.dim,arg.split);
[y,yobs,yvar] = formatcells(y,arg.dim,arg.split);
[z1,zobs1] = formatcells(resp1,arg.dim,arg.split);
[z2,zobs2] = formatcells(resp2,arg.dim,arg.split);

% Check equal number of observations
if ~isequal(xobs,yobs,zobs1,zobs2)
    error(['STIM and RESP arguments must have the same number of '...
        'observations.'])
end

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*Dir);
tmax = ceil(tmax/1e3*fs*Dir);
lags = tmin:tmax;
arg.window = round(arg.window*fs);

% Compute sampling interval
delta = 1/fs;

% Get dimensions
xvar = unique(xvar);
yvar = unique(yvar);
nreg = numel(lambda);
nfold = numel(x);
switch arg.type
    case 'multi'
        nvar = xvar*numel(lags)+1;
        nlag = 1;
    case 'single'
        nvar = xvar+1;
        nlag = numel(lags);
end

% Truncate output
if ~arg.zeropad
    [y,yobs] = truncate(y,tmin,tmax,yobs);
    if Dir == 1
        z1 = truncate(z1,tmin,tmax,zobs1);
        z2 = truncate(z2,tmin,tmax,zobs2);
    end
end
if arg.window
    nwin = sum(floor(yobs/arg.window));
else
    nwin = nfold;
end

% Verbose mode
if arg.verbose
    v = verbosemode([],[],nfold);
end

% Compute unisensory covariance matrices
if Dir == 1
    if arg.fast
        [Cxx,Cxy1,Cxy2,folds] = mlscovmat(x,z1,z2,lags,arg.type,...
            arg.zeropad,arg.verbose);
    else
        [Cxx,Cxy1,Cxy2] = mlscovmat(x,z1,z2,lags,arg.type,...
            arg.zeropad,arg.verbose);
    end
elseif Dir == -1
    if arg.fast
        [Cxx1,Cxy1,folds] = olscovmat(z1,y,lags,arg.type,...
            arg.zeropad,arg.verbose);
        [Cxx2,Cxy2,folds2] = olscovmat(z2,y,lags,arg.type,...
            arg.zeropad,arg.verbose);
    else
        [Cxx1,Cxy1] = olscovmat(z1,y,lags,arg.type,...
            arg.zeropad,arg.verbose);
        [Cxx2,Cxy2] = olscovmat(z2,y,lags,arg.type,...
            arg.zeropad,arg.verbose);
    end
end

% Verbose mode
if arg.verbose
    v = verbosemode(v,0,nfold);
end

% Sum covariances for additive model
if Dir == 1
    Cxx = Cxx + Cxx;
    Cxy = Cxy1 + Cxy2;
    if arg.fast
        for i = 1:nfold
            folds.Cxx{i} = folds.Cxx{i} + folds.Cxx{i};
            folds.Cxy{i} = folds.Cxy{i} + folds.Cxz{i};
        end
    end
elseif Dir == -1
    Cxx = Cxx1 + Cxx2;
    Cxy = Cxy1 + Cxy2;
    if arg.fast
        for i = 1:nfold
            folds.Cxx{i} = folds.Cxx{i} + folds2.Cxx{i};
            folds.Cxy{i} = folds.Cxy{i} + folds2.Cxy{i};
        end
    end
end

% Set up sparse regularization matrix
M = regmat(nvar,arg.method)/delta;

% Initialize variables
r = zeros(nwin,nreg,yvar,nlag);
err = zeros(nwin,nreg,yvar,nlag);
ii = 0;

% Leave-one-out cross-validation
for i = 1:nfold
    
    if arg.window
        ii = ii(end)+1:ii(end)+floor(yobs(i)/arg.window);
    else
        ii = i;
    end
    
    if arg.fast % fast method
        
        % Multisensory validation set
        if Dir == 1
            xlag = folds.xlag{i};
        elseif Dir == -1
            xlag = lagGen(x{i},lags,arg.zeropad,1);
        end
        
        % Unisensory training set
        Cxxi = Cxx - folds.Cxx{i};
        Cxyi = Cxy - folds.Cxy{i};
        
    else % memory-efficient method
        
        % Multisensory validation set
        xlag = lagGen(x{i},lags,arg.zeropad,1);
        
        if Dir == 1
            
            % Unisensory training set
            Cxxi = Cxx - xlag'*xlag*2;
            Cxyi = Cxy - xlag'*y{i}*2;
            
        elseif Dir == -1
            
            % Unisensory left-out set
            z1lag = lagGen(z1{i},lags,arg.zeropad,1);
            z2lag = lagGen(z2{i},lags,arg.zeropad,1);
            
            % Unisensory training set
            Cxxi = Cxx - (z1lag'*z1lag + z2lag'*z2lag);
            Cxyi = Cxy - (z1lag'*y{i} + z2lag'*y{i});
            
        end
        
    end
    
    for j = 1:nreg
        
        switch arg.type
            
            case 'multi'
                
                % Fit linear model
                w = (Cxxi + lambda(j)*M)\Cxyi*2;
                
                % Predict multisensory output
                pred = xlag*w;
                
                % Evaluate performance
                [r(ii,j,:),err(ii,j,:)] = mTRFevaluate(y{i},pred,...
                    'corr',arg.corr,'error',arg.error,...
                    'window',arg.window);
                
            case 'single'
                
                for k = 1:nlag
                    
                    % Index lag
                    idx = [1,xvar*(k-1)+2:xvar*k+1];
                    
                    % Fit linear model
                    w = (Cxxi(:,:,k) + lambda(j)*M)\Cxyi(:,:,k)*2;
                    
                    % Predict multisensory output
                    pred = xlag(:,idx)*w;
                    
                    % Evaluate performance
                    [r(ii,j,:,k),err(ii,j,:,k)] = mTRFevaluate(y{i},pred,...
                        'corr',arg.corr,'error',arg.error,...
                        'window',arg.window);
                    
                end
                
        end
        
    end
    
    % Verbose mode
    if arg.verbose
        v = verbosemode(v,i,nfold);
    end
    
end

% Format output arguments
stats = struct('r',r,'err',err);
if nargout > 1
    t = lags/fs*1e3;
end

% Verbose mode
if arg.verbose
    verbosemode(v,i+1,nfold,stats);
end

function v = verbosemode(v,fold,nfold,stats)
%VERBOSEMODE  Execute verbose mode.
%   V = VERBOSEMODE(V,FOLD,NFOLD,STATS) prints details about the progress
%   of the main function to the screen.

if isempty(fold)
    v = struct('msg',[],'h',[],'tocs',0);
    fprintf('\nTrain on %d folds, validate on 1 fold\n',nfold-1)
elseif fold == 0
    fprintf('Training/validating model\n')
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
    fprintf('val_correlation: %.4f - val_error: %.4f\n',rmax,emax)
end

function validateparamin(fs,Dir,tmin,tmax,lambda)
%VALIDATEPARAMIN  Validate input parameters.
%   VALIDATEPARAMIN(FS,DIR,TMIN,TMAX,LAMBDA) validates the input parameters
%   of the main function.

if ~isnumeric(fs) || ~isscalar(fs) || fs <= 0
    error('FS argument must be a positive numeric scalar.')
elseif Dir ~= 1 && Dir ~= -1
    error('DIR argument must have a value of 1 or -1.')
elseif ~isnumeric([tmin,tmax]) || ~isscalar(tmin) || ~isscalar(tmax)
    error('TMIN and TMAX arguments must be numeric scalars.')
elseif tmin > tmax
    error('The value of TMIN must be less than that of TMAX.')
elseif ~isnumeric(lambda) || any(lambda < 0)
    error('LAMBDA argument must be positive numeric values.')
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

% Regularization method
regOptions = {'ridge','Tikhonov','ols'};
validFcn = @(x) any(validatestring(x,regOptions));
addParameter(p,'method','ridge',validFcn);

% Model type
lagOptions = {'multi','single'};
validFcn = @(x) any(validatestring(x,lagOptions));
addParameter(p,'type','multi',validFcn);

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
addParameter(p,'fast',true,validFcn); % fast CV method
addParameter(p,'verbose',true,validFcn); % verbose mode

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially matched strings
arg.method = validatestring(arg.method,regOptions);
arg.type = validatestring(arg.type,lagOptions);
arg.corr = validatestring(arg.corr,corrOptions);
arg.error = validatestring(arg.error,errOptions);