function [stats,stats1,stats2,t] = mTRFattncrossval(stim1,stim2,resp,fs,Dir,tmin,tmax,lambda,varargin)
%MTRFATTNCROSSVAL  Cross-validation for attention decoder optimization.
%   STATS = MTRFATTNCROSSVAL(STIM1,STIM2,RESP,FS,DIR,TMIN,TMAX,LAMBDA)
%   cross validates a forward encoding model (stimulus to neural response)
%   or a backward decoding model (neural response to stimulus) over
%   multiple trials of data for optimizing an attention decoder. Models are
%   trained on the attended stimuli STIM1, and validated on both STIM1 and
%   the unattended stimuli STIM2 as per O'Sullivan et al. (2015). Pass in 1
%   for DIR to validate a forward model, or -1 to validate a backward
%   model. STIM1, STIM2 and RESP are cell arrays containing corresponding
%   trials of continuous data. FS is a scalar specifying the sample rate in
%   Hertz, and TMIN and TMAX are scalars specifying the minimum and maximum
%   time lags in milliseconds. For backward models, the time lags are
%   automatically reversed. LAMBDA is a vector of regularization values to
%   be validated and controls overfitting.
%
%   MTRFATTNCROSSVAL returns the cross-validation statistics of the
%   attention decoder in a structure with the following fields:
%       'acc'       -- decoding accuracy based on the proportion of
%                      observations where the correlation for STIM1 was
%                      greater than that of STIM2 (nlambda-by-yvar)
%       'd'         -- attention modulation index based on d', where the
%                      correlation for STIM1 is considered signal and that
%                      of STIM2 is considered noise (nlambda-by-yvar)
%
%   MTRFATTNCROSSVAL performs a leave-one-out cross-validation over all
%   trials. To achieve a k-fold cross-validation, arrange STIM1, STIM2 and
%   RESP in k-by-1 cell arrays. The number of folds can also be increased
%   by an integer factor using the 'split' parameter (see below).
%
%   If STIM1, STIM2 or RESP contain matrices, it is assumed that the rows
%   correspond to observations and the columns to variables, unless
%   otherwise stated via the 'dim' parameter (see below). If they contain
%   vectors, it is assumed that the first non-singleton dimension
%   corresponds to observations. Each trial of STIM1, STIM2 and RESP must
%   have the same number of observations.
%
%   [STATS,STATS1,STATS2] = MTRFATTNCROSSVAL(...) returns the cross-
%   validation statistics for STIM1 and STIM2 respectively in structures
%   with the following fields:
%       'r'         -- correlation coefficient based on Pearson's linear
%                      correlation coefficient (nfold-by-nlambda-by-yvar)
%       'err'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar)
%
%   [STATS,STATS1,STATS2,T] = MTRFATTNCROSSVAL(...) returns a vector
%   containing the time lags used in milliseconds. These data are useful
%   for interpreting the results of single-lag models.
%
%   [...] = MTRFATTNCROSSVAL(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
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
%   Each iteration of MTRFATTNCROSSVAL partitions the N trials or segments
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
%   See also CROSSVAL, DPRIME, MTRFPARTITION, MTRFTRAIN, MTRFPREDICT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] O'Sullivan JA, Power AJ, Mesgarani N, Rajaram S, Foxe JJ, Shinn-
%          Cunningham BG, Slaney M, Shamma SA, Lalor EC (2015) Attentional
%          Selection in a Cocktail Party Environment Can Be Decoded from
%          Single-Trial EEG. Cereb Cortex 25(7):1697-1706.
%      [2] Alickovic E, Lunner T, Gustafsson F, Ljung L (2019) A Tutorial
%          on Auditory Attention Identification Methods. Front Neurosci
%          13:153.

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
    x = stim1; y = resp;
elseif Dir == -1
    x = resp; y = stim1;
    [tmin,tmax] = deal(tmax,tmin);
end

% Format data in cell arrays
[x,xobs,xvar] = formatcells(x,arg.dim,arg.split);
[y,yobs,yvar] = formatcells(y,arg.dim,arg.split);
[z,zobs] = formatcells(stim2,arg.dim,arg.split);

% Check equal number of observations
if ~isequal(xobs,yobs,zobs)
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
    if Dir == -1
        z = truncate(z,tmin,tmax,zobs);
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

% Compute covariance matrices
if arg.fast
    [Cxx,Cxy,folds] = olscovmat(x,y,lags,arg.type,arg.zeropad,arg.verbose);
else
    [Cxx,Cxy] = olscovmat(x,y,lags,arg.type,arg.zeropad,arg.verbose);
end

% Verbose mode
if arg.verbose
    v = verbosemode(v,0,nfold);
end

% Set up sparse regularization matrix
M = regmat(nvar,arg.method)/delta;

% Initialize variables
r1 = zeros(nwin,nreg,yvar,nlag);
err1 = zeros(nwin,nreg,yvar,nlag);
r2 = zeros(nwin,nreg,yvar,nlag);
err2 = zeros(nwin,nreg,yvar,nlag);
ii = 0;

% Leave-one-out cross-validation
for i = 1:nfold
    
    if arg.window
        ii = ii(end)+1:ii(end)+floor(yobs(i)/arg.window);
    else
        ii = i;
    end
    
    if arg.fast % fast method
        
        % Validation set
        xlag = folds.xlag{i};
        
        % Training set
        Cxxi = Cxx - folds.Cxx{i};
        Cxyi = Cxy - folds.Cxy{i};
        
    else % memory-efficient method
        
        % Validation set
        xlag = lagGen(x{i},lags,arg.zeropad);
        
        % Training set
        Cxxi = Cxx - xlag'*xlag;
        Cxyi = Cxy - xlag'*y{i};
        
    end
    
    % Unattended validation set
    if Dir == 1
        zlag = lagGen(z{i},lags,arg.zeropad,1);
    end
    
    for j = 1:nreg
        
        switch arg.type
            
            case 'multi'
                
                % ---Attended Stimulus---
                
                % Fit linear model
                w = (Cxxi + lambda(j)*M)\Cxyi;
                
                % Predict output
                pred = xlag*w;
                
                % Evaluate performance
                [r1(ii,j,:),err1(ii,j,:)] = mTRFevaluate(y{i},pred,...
                    'corr',arg.corr,'error',arg.error,...
                    'window',arg.window);
                
                % ---Unattended Stimulus---
                
                if Dir == 1
                    
                    % Predict output
                    pred = zlag*w;
                    
                    % Evaluate performance
                    [r2(ii,j,:),err2(ii,j,:)] = mTRFevaluate(y{i},pred,...
                        'corr',arg.corr,'error',arg.error,...
                        'window',arg.window);
                    
                elseif Dir == -1
                    
                    % Evaluate performance
                    [r2(ii,j,:),err2(ii,j,:)] = mTRFevaluate(z{i},pred,...
                        'corr',arg.corr,'error',arg.error,...
                        'window',arg.window);
                    
                end
                
            case 'single'
                
                for k = 1:nlag
                    
                    % Index lag
                    idx = [1,xvar*(k-1)+2:xvar*k+1];
                    
                    % ---Attended Stimulus---
                    
                    % Fit linear model
                    w = (Cxxi(:,:,k) + lambda(j)*M)\Cxyi(:,:,k);
                    
                    % Predict output
                    pred = xlag(:,idx)*w;
                    
                    % Evaluate performance
                    [r1(ii,j,:,k),err1(ii,j,:,k)] = ...
                        mTRFevaluate(y{i},pred,...
                        'corr',arg.corr,'error',arg.error,...
                        'window',arg.window);
                    
                    % ---Unattended Stimulus---
                    
                    if Dir == 1
                        
                        % Predict output
                        pred = zlag(:,idx)*w;
                        
                        % Evaluate performance
                        [r2(ii,j,:,k),err2(ii,j,:,k)] = ...
                            mTRFevaluate(y{i},pred,...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                        
                    elseif Dir == -1
                        
                        % Evaluate performance
                        [r2(ii,j,:,k),err2(ii,j,:,k)] = ...
                            mTRFevaluate(z{i},pred,...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                        
                    end
                    
                end
                
        end
        
    end
    
    % Verbose mode
    if arg.verbose
        v = verbosemode(v,i,nfold);
    end
    
end

% Compute decoding accuracy and d'
[acc,d] = mTRFattnevaluate(r1,r2);

% Format output arguments
stats = struct('acc',acc,'d',d);
stats1 = struct('r',r1,'err',err1);
stats2 = struct('r',r2,'err',err2);
if nargout > 3
    t = lags/fs*1e3;
end

% Verbose mode
if arg.verbose
    verbosemode(v,i+1,nfold,stats,stats1,stats2);
end

function v = verbosemode(v,fold,nfold,stats,stats1,stats2)
%VERBOSEMODE  Execute verbose mode.
%   V = VERBOSEMODE(V,FOLD,NFOLD,STATS,STATS1,STATS2) prints details about
%   the progress of the main function to the screen.

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
    rmax = mean(stats1.r,1); rmax = max(rmax(:));
    emax = mean(stats1.err,1); emax = max(emax(:));
    fprintf('STIM1 - val_correlation: %.4f - val_error: %.4f\n',rmax,emax)
    rmax = mean(stats2.r,1); rmax = max(rmax(:));
    emax = mean(stats2.err,1); emax = max(emax(:));
    fprintf('STIM2 - val_correlation: %.4f - val_error: %.4f\n',rmax,emax)
    amax = max(stats.acc(:)); dmax = max(stats.d(:));
    fprintf('S1vS2 - val_accuracy: %.4f - val_d'': %.4f\n',amax,dmax)
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