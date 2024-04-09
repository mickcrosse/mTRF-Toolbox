function model = mTRFtrain(stim,resp,fs,Dir,tmin,tmax,lambda,varargin)
%MTRFTRAIN  Train a linear encoding/decoding model.
%   MODEL = MTRFTRAIN(STIM,RESP,FS,DIR,TMIN,TMAX,LAMBDA) trains a forward
%   encoding model (stimulus to neural response) or a backward decoding
%   model (neural response to stimulus) using time-lagged input features as
%   per Crosse et al. (2016). Pass in 1 for DIR to fit a forward model, or
%   -1 to fit a backward model. STIM and RESP are matrices or cell arrays
%   containing corresponding trials of continuous training data. FS is a
%   scalar specifying the sample rate in Hertz, and TMIN and TMAX are
%   scalars specifying the minimum and maximum time lags in milliseconds.
%   For backward models, the time lags are automatically reversed. LAMBDA
%   is a scalar specifying the regularization parameter for controlling
%   overfitting.
%
%   MTRFTRAIN returns the model in a structure with the following fields:
%       'w'         -- normalized model weights (xvar-by-nlag-by-yvar)
%       'b'         -- normalized bias term (1-by-nlag-by-yvar)
%       't'         -- time lags (ms)
%       'fs'        -- sample rate (Hz)
%       'Dir'       -- direction of causality (forward=1, backward=-1)
%       'type'      -- type of model (multi-lag, single-lag)
%
%   MTRFTRAIN normalizes the model weights and regularization matrix by the
%   sampling interval (1/FS) such that their magnitude and smoothness is
%   invariant to sample rate FS (Lalor et al., 2006).
%
%   If STIM or RESP are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   STIM and RESP must have the same number of observations.
%
%   If STIM and RESP are cell arrays containing multiple trials, the
%   covariance matrices of each trial are summed to produce one model. STIM
%   and RESP must contain the same number of trials.
%
%   [...] = MTRFTRAIN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%       'split'     A scalar specifying the number of segments in which to
%                   split each trial of data when computing the covariance
%                   matrices. This is useful for reducing memory usage on
%                   large datasets. By default, the entire trial is used.
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them. Pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%       'verbose'   A numeric or logical specifying whether to execute in
%                   verbose mode: pass in 1 for verbose mode (default), or
%                   0 for non-verbose mode.
%
%   See also RIDGE, REGRESS, MTRFTRANSFORM, MTRFCROSSVAL.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Lalor EC, Pearlmutter BA, Reilly RB, McDarby G, Foxe JJ (2006)
%          The VESPA: a method for the rapid estimation of a visual evoked
%          potential. NeuroImage 32:1549-1561.

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

% Check equal number of observations
if ~isequal(xobs,yobs)
    error(['STIM and RESP arguments must have the same number of '...
        'observations.'])
end

% Verbose mode
if arg.verbose
    v = verbosemode(1,sum(xobs));
end

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*Dir);
tmax = ceil(tmax/1e3*fs*Dir);
lags = tmin:tmax;

% Compute sampling interval
delta = 1/fs;

% Get dimensions
nlag = numel(lags);
xvar = unique(xvar);
yvar = unique(yvar);
switch arg.type
    case 'multi'
        nvar = xvar*nlag+1;
    case 'single'
        nvar = xvar+1;
end

% Truncate output
if ~arg.zeropad
    y = truncate(y,tmin,tmax,yobs);
end

% Compute covariance matrices
[Cxx,Cxy] = olscovmat(x,y,lags,arg.type,arg.zeropad,arg.verbose);

% Verbose mode
if arg.verbose
    v = verbosemode(v);
end

% Set up sparse regularization matrix
M = regmat(nvar,arg.method)*lambda/delta;

% Fit linear model
switch arg.type
    case 'multi'
        w = (Cxx + M)\Cxy/delta;
    case 'single'
        w = zeros(nvar,nlag,yvar);
        for i = 1:nlag
            w(:,i,:) = (Cxx(:,:,i) + M)\Cxy(:,:,i)/delta;
        end
end

% Format output arguments
model = struct('w',reshape(w(2:end,:,:),[xvar,nlag,yvar]),'b',w(1,:,:),...
    't',lags/fs*1e3,'fs',fs,'Dir',Dir,'type',arg.type);

% Verbose mode
if arg.verbose
    verbosemode(v,[],model);
end

function v = verbosemode(v,nobs,model)
%VERBOSEMODE  Execute verbose mode.
%   V = VERBOSEMODE(V,NOBS,MODEL) prints details about the progress of the
%   main function to the screen.

if v == 1
    fprintf('\nTrain on %d samples\n',nobs)
elseif v == 2
    fprintf('Training model'); tic
elseif v == 3
    fprintf(' - %.3fs\n',toc)
    modelsummary(model)
end

v = v+1;

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
elseif ~isnumeric(lambda) || ~isscalar(lambda) || lambda < 0
    error('LAMBDA argument must be a positive numeric scalar.')
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

% Split data
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'split',1,validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'zeropad',true,validFcn); % zero-pad design matrix
addParameter(p,'verbose',true,validFcn); % verbose mode

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially-matched strings
arg.method = validatestring(arg.method,regOptions);
arg.type = validatestring(arg.type,lagOptions);