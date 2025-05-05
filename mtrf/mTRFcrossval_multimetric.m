function [stats,t] = mTRFcrossval_multimetric(stim,resp,fs,Dir,tmin,tmax,lambda,varargin)

% MTRFCROSSVAL  Leave-one-out cross-validation.
%   STATS = MTRFCROSSVAL_MULTIMETRIC(STIM,RESP,GROUNDTRUTH,FS,DIR,TMIN,TMAX,LAMBDA,EVALTYPE,WINST,WINEND,WODIM) cross validates
%   a forward encoding model (stimulus to neural response) or a backward
%   decoding model (neural response to stimulus) over multiple trials of
%   data as per Crosse et al. (2016). Pass in 1 for DIR to validate a
%   forward model, or -1 to validate a backward model. STIM and RESP are
%   cell arrays containing corresponding trials of continuous data.
%	GROUNDTRUTH is the average responses across subjects for every trial.	
% 	FS is a scalar specifying the sample rate in Hertz, and TMIN and TMAX are
%   scalars specifying the minimum and maximum time lags in milliseconds.
%   For backward models, the time lags are automatically reversed. LAMBDA
%   is a vector of regularization values to be validated and controls
%   overfitting.
%
%   In addition to the standard mTRFcrossval, this function can return
%   additional evaluation metrics. The first metric, named GT, fits the TRF
%   models as usual, but evaluates their performance on a "cleaner" neural
%   signal, which we call ground truth. This signal might be the actual
%   ground truth, which is available when working with simulated data, or
%   an estimate of the signal buried behind the EEG noise e.g., average
%   across data from many individuals, MCCA components (see Alain
%   de Cheveigne's work). The second metric is called TC, which stands for
%   time-constrained. As for the other metric, TRF models are fit as usual,
%   while only the evaluation changes. Specifically, TC calculates the
%   prediction correlations only on certain samples of interest. This is
%   useful, for example, when deadling with event onsets such as words,
%   where the effect of interest might only be between 300 and 500ms after the
%   onset of each word. Calculating correlations on all the samples would
%   dilute the effect. As such, this function enables the user to specify
%   which dimension indicates the onsets of interest as well as the window
%   of interest according to those onsets.
%
%   MTRFCROSSVAL returns the cross-validation statistics in a structure
%   with some or all following fields, depending on the parameters:
%       'r'         -- correlation coefficient based on Pearson's linear
%                      correlation coefficient (nfold-by-nlambda-by-yvar) for the original mTRFcrossval metric
%       'err'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar)for the original mTRFcrossval metric
%
%       'r_TC'         -- correlation coefficient based on Pearson's linear
%                      correlation coefficient (nfold-by-nlambda-by-yvar) for the TC metric 
%       'err_TC'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar)for the TC metric metric
%
%       'r_GT_TC'         -- correlation coefficient based on Pearson's linear
%                      correlation coefficient (nfold-by-nlambda-by-yvar)
%                      for the the GT and TC metrics
%       'err_GT_TC'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar)for he GT and TC metrics
%
%       'r_GT'         -- correlation coefficient based on Pearson's linear 
%                      correlation coefficient (nfold-by-nlambda-by-yvar) for the GT metric
%       'err_GT'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar) for the GT metric
%
%   MTRFCROSSVAL performs a leave-one-out cross-validation over all trials.
%   To achieve a k-fold cross-validation, arrange STIM and RESP in k-by-1
%   cell arrays. The number of folds can also be increased by an integer
%   factor using the 'split' parameter (see below).
%
%   If STIM or RESP contain matrices, it is assumed that the rows
%   correspond to observations and the columns to variables, unless
%   otherwise stated via the 'dim' parameter (see below). If they contain
%   vectors, it is assumed that the first non-singleton dimension
%   corresponds to observations. Each trial of STIM and RESP must have the
%   same number of observations.
%
%   [STATS,T] = MTRFCROSSVAL(...) returns a vector containing the time lags
%   used in milliseconds. These data are useful for interpreting the
%   results of single-lag models.
%
%   [...] = MTRFCROSSVAL(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%
%       'respGroundTruth'  a version of 'resp' containing the ground truth
%                   signal buried behind the noise. Such a signal can be
%                   available when working with simulated neural data. When
%                   working with real data, the ground-truth can be
%                   roughly estimated (with numerous caveats) with
%                   methodologies such as MCCA (see Alain de Cheveigne's
%                   work) or by simply averaging the neural signals across
%                   many individuals (again, with strong assumptions).
%
%       'onsetVecDim'   scalar indicating the dimension of stim.data 
%                   matrix corresponding to the onset vector (e.g., word 
%                   onset, word surprisal) in your multivariate 
%       'onsetVecWinSt' and 'onsetVecWinEnd' are scalars indicating the
%                   beginning and the end of the window for the 'TC" and
%                   'GT_TC' metric in milliseconds 
%       
%
%       'verbose'   A numeric or logical specifying whether to execute in
%                   verbose mode: pass in 1 for verbose mode (default), or
%                   0 for non-verbose mode.
%
%   Notes:
%   Each iteration of MTRFCROSSVAL partitions the N trials or segments of
%   data into two subsets, fitting a model to N-1 trials (training set) and
%   validating it on the left-out trial (validation set). Performance on
%   the validation set can be used to optimize hyperparameters (e.g.,
%   LAMBDA, TMAX). Once completed, it is recommended to evaluate model
%   performance on separate held-out data using the mTRFpredict function.
%
%   Discontinuous trials of data should not be concatenated prior to cross-
%   validation, as this will introduce artifacts in places where the
%   temporal integration window crosses over trial boundaries. Each trial
%   of continuous data should be input as a separate cell.
%
%   Example: [stats,t] = mTRFcrossval_multimetric(stimFeature.data,eeg.data, ...
%             eeg.fs,dirTRF,tmin,tmax,lambdas,'verbose',0,'onsetVecDim',1, ...
%             'onsetVecWinSt',300,'onsetVecWinEnd',500,'respGroundTruth',eeg.data);
%   
%   See also CROSSVAL, MTRFPARTITION, MTRFTRAIN, MTRFPREDICT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Alickovic E, Lunner T, Gustafsson F, Ljung L (2019) A Tutorial
%          on Auditory Attention Identification Methods. Front Neurosci
%          13:153.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Edmund Lalor <edmundlalor@gmail.com>
%            Nate Zuk <zukn@tcd.ie>
%            Emily Ip <ipem@tcd.ie>
%            Amirhossein Chalehchaleh <chalehca@tcd.ie>
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

if iscell(arg.respGroundTruth)
    [respGroundTruth,~,~] = formatcells(arg.respGroundTruth,arg.dim,arg.split);
end

% Check equal number of observations
if ~isequal(xobs,yobs)
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
r = zeros(nwin,nreg,yvar,nlag);
err = zeros(nwin,nreg,yvar,nlag);

r_GT = zeros(nwin,nreg,yvar,nlag);
err_GT = zeros(nwin,nreg,yvar,nlag);

r_GT_TC = zeros(nwin,nreg,yvar,nlag);
err_GT_TC = zeros(nwin,nreg,yvar,nlag);

r_TC = zeros(nwin,nreg,yvar,nlag);
err_TC = zeros(nwin,nreg,yvar,nlag);

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
        xlag = lagGen(x{i},lags,arg.zeropad,1);

        % Training set
        Cxxi = Cxx - xlag'*xlag;
        Cxyi = Cxy - xlag'*y{i};

    end

    if arg.onsetVecDim > 0 % if TC metric in use
        % Define interval for computing TC correlation
        st = round(arg.onsetVecWinSt/1000*fs);
        fin = round(arg.onsetVecWinEnd/1000*fs);
        interval = st:fin;

        % Find index of word onsets
        ind1 = find(x{i}(:,arg.onsetVecDim)~=0);

        if isempty(ind1) % empty trial or no onset in trial
            continue
        end

        idx=[];
        for iin=interval
            idx=cat(2,idx,ind1+iin);
        end
        idx=reshape(transpose(idx),[],1);

        idx=unique(idx); % unique indices only

        if any(idx>length(x{i}))
            % window size larger than stim size (last onset)
            [row, col]=find(idx>length(x{i}));

            idx(row)=length(x{i});
            idx=unique(idx);
        end
    end

    % Evaluating performance
    for j = 1:nreg

        switch arg.type

            case 'multi'

                % Fit linear model
                w = (Cxxi + lambda(j)*M)\Cxyi;

                % Predict output
                pred = xlag*w;

                % Standard metric
                [r(ii,j,:),err(ii,j,:)] = mTRFevaluate(y{i},pred,...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                % TC
                if arg.onsetVecDim > 0
                    [r_TC(ii,j,:),err_TC(ii,j,:)] = mTRFevaluate(y{i}(idx,:),pred(idx,:),...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                end
                
                % GT
                if iscell(arg.respGroundTruth)
                    [r_GT(ii,j,:),err_GT(ii,j,:)] = mTRFevaluate(respGroundTruth{i},pred,...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                end
                
                % GT_TC
                if arg.onsetVecDim > 0 && iscell(arg.respGroundTruth)
                    [r_GT_TC(ii,j,:),err_GT_TC(ii,j,:)] = mTRFevaluate(respGroundTruth{i}(idx,:),pred(idx,:),...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                end

            case 'single'

                for k = 1:nlag

                    % Index lag
                    idx = [1,xvar*(k-1)+2:xvar*k+1];

                    % Fit linear model
                    w = (Cxxi(:,:,k) + lambda(j)*M)\Cxyi(:,:,k);

                    % Predict output
                    pred = xlag(:,idx)*w;

                    % Evaluate performance
                    % Standard metric
                    [r(ii,j,:,k),err(ii,j,:,k)] = mTRFevaluate(y{i},pred,...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                    
                    % TC
                    if arg.onsetVecDim > 0
                        [r_GT_TC(ii,j,:),err_GT_TC(ii,j,:)] = mTRFevaluate(y{i}(idx,:),pred(idx,:),...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                    end
                    
                    % GT
                    if iscell(arg.respGroundTruth)
                        [r(ii,j,:,k),err(ii,j,:,k)] = mTRFevaluate(respGroundTruth{i},pred,...
                            'corr',arg.corr,'error',arg.error,...
                            'window',arg.window);
                    end
                    
                    % GT_TC
                    if arg.onsetVecDim > 0 && iscell(arg.groundTruth)
                        [r_GT_TC(ii,j,:),err_GT_TC(ii,j,:)] = mTRFevaluate(respGroundTruth{i}(idx,:),pred(idx,:),...
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

    if iscell(arg.respGroundTruth)
        predAll{i} = pred;
    end

end

% Format output arguments
% save all eval metric outputs
stats = struct('r',r,'err',err, ...             % OG
               'r_TC',r_TC,'err_TC',err_TC, ...  % TC
               'r_GT_TC',r_GT_TC,'err_GT_TC',err_GT_TC,  ...
               'r_GT', r_GT, 'err_GT', err_GT); % GT

stats = struct('r',r,'err',err);
if arg.onsetVecDim > 0 % TC
    stats.r_TC = r_TC;
    stats.err_TC = err_TC;
end

if iscell(arg.respGroundTruth) % GT
    stats.r_GT = r_GT;
    stats.err_GT = err_GT;
end

if arg.onsetVecDim > 0 && iscell(arg.respGroundTruth) % GT_TC
    stats.r_GT_TC = r_GT_TC;
    stats.err_GT_TC = err_GT_TC;
end
                        
if iscell(arg.respGroundTruth)
    stats.predAll = predAll;
end

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

% respGroundTruth
errorMsg = 'It must be a cell.';
validFcn = @(x) assert(iscell(x),errorMsg);
addParameter(p,'respGroundTruth',0,validFcn);

% onsetVecDim
errorMsg = 'It must be an integer greater than zero.';
validFcn = @(x) assert((fix(x)-x)==0 && x>=1,errorMsg);
addParameter(p,'onsetVecDim',0,validFcn);

% onsetVecWinSt
errorMsg = 'It must be a numerical value.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'onsetVecWinSt',0,validFcn);

% onsetVecWinEnd
errorMsg = 'It must be a numerical value.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'onsetVecWinEnd',0,validFcn);

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially matched strings
arg.method = validatestring(arg.method,regOptions);
arg.type = validatestring(arg.type,lagOptions);
arg.corr = validatestring(arg.corr,corrOptions);
arg.error = validatestring(arg.error,errOptions);