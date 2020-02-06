function [r,p,rmse,t] = mTRFcrossval(stim,resp,fs,dir,tmin,tmax,lambda,varargin)
%MTRFCROSSVAL  mTRF-Toolbox cross-validation.
%   [R,P,RMSE] = MTRFCROSSVAL(STIM,RESP,FS,DIR,TMIN,TMAX,LAMBDA) cross
%   validates a forward encoding model (stimulus to neural response) or a
%   backward decoding model (neural response to stimulus) using time-lagged
%   input features. Pass in 1 for DIR to validate a forward model, or -1 to
%   validate a backward model. STIM and RESP are cell arrays containing
%   corresponding trials of continuous data over which to cross-validate.
%   FS is a scalar specifying the sample rate in Hertz, and TMIN and TMAX
%   are scalars specifying the minimum and maximum time lags in
%   milliseconds. For backward models, MTRFCROSSVAL automatically reverses
%   the time lags. LAMBDA is a scalar or vector of regularization values
%   to be validated and controls overfitting.
%
%   As a measure of performance, MTRFCROSSVAL returns matrices containing
%   the correlation coefficients R between the predicted and observed
%   variables, the the probabilities of the correlation coefficients P and
%   the root-mean-square errors RMSE between the predicted and observed
%   variables. The size of each matrix is ntrial-by-nlambda-by-yvar. For
%   single-lag models, a fourth dimension is added to specify the lag.
%
%   MTRFCROSSVAL performs a leave-one-out cross-validation across trials.
%   To achieve a k-fold cross-validation, arrange STIM and RESP in k-by-1
%   cell arrays. The number of folds can also be increased by an integer
%   factor using the 'split' parameter (see below). It is not recommended
%   to concatenate discontinuous data segments or to use cross-validation
%   as a way to test model performance. Models should be tested on separate
%   data after cross-validation using the mTRFpredict function.
%
%   If STIM or RESP contain matrices, it is assumed that the rows
%   correspond to observations and the columns to variables, unless
%   otherwise stated via the 'dim' parameter (see below). If they contain
%   vectors, it is assumed that the first non-singleton dimension
%   corresponds to observations. Each trial of STIM and RESP must have the
%   same number of observations.
%
%   [R,P,RMSE,T] = MTRFCROSSVAL(...) returns a vector containing the time
%   lags used in milliseconds. This is useful for plotting single-lag model
%   results.
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
%       'type'      A string specifying type of model to fit:
%                       'multi'     use all lags simultaneously to fit a
%                                   multi-lag model (default)
%                       'single'    use each lag individually to fit
%                                   separate single-lag models
%       'split'     A scalar specifying the number of segments in which to
%                   split each trial of data when computing the covariance
%                   matrices. This is useful for reducing memory usage on
%                   large datasets. By default, the entire trial is used.
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%       'fast'      A numeric or logical specifying whether to use the fast
%                   cross-validation method (requires more memory) or the
%                   slower method (requires less memory): pass in 1 to use
%                   the fast method (default) or 0 to use the slower
%                   method. Note, both methods are numerically equivalent.
%
%   See mTRFdemos for examples of use.
%
%   See also MTRFTRAIN, MTRFPREDICT, MTRFTRANSFORM, MTRFMULTICROSSVAL.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2015) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Alickovic E, Lunner T, Gustafsson F, Ljung L (2019) A Tutorial
%          on Auditory Attention Identification Methods. Front Neurosci
%          13:153.

%   Authors: Mick Crosse, Giovanni Di Liberto, Nate Zuk, Edmund Lalor
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Apr 2014; Last revision: 05-Feb-2020

% Parse input arguments
arg = parsevarargin(varargin);

% Validate parameter values
if ~isnumeric(fs) || ~isscalar(fs) || fs <= 0
    error('FS argument must be a positive numeric scalar.')
elseif ~isnumeric([tmin,tmax]) || ~isscalar(tmin) || ~isscalar(tmax)
    error('TMIN and TMAX arguments must be numeric scalars.')
elseif tmin > tmax
    error('The value of TMIN must be less than that of TMAX.')
elseif ~isnumeric(lambda) || any(lambda < 0)
    error('LAMBDA argument must be positive numeric values.')
end

% Define X and Y variables
if dir == 1
    x = stim; y = resp;
elseif dir == -1
    x = resp; y = stim;
    [tmin,tmax] = deal(tmax,tmin);
else
    error('DIR argument must have a value of 1 or -1.')
end

% Format data in cells column-wise
[x,xobs,xvar] = formatcells(x,arg.dim);
[y,yobs,yvar] = formatcells(y,arg.dim);

% Check equal number of observations
if ~isequal(xobs,yobs)
    error(['STIM and RESP arguments must have the same number of '...
        'observations.'])
end

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*dir);
tmax = ceil(tmax/1e3*fs*dir);
lags = tmin:tmax;

% Compute sampling interval
delta = 1/fs;

% Get dimensions
nlag = numel(lags);
xvar = unique(xvar);
yvar = unique(yvar);
nvar = xvar*nlag+1;
switch arg.type
    case 'multi'
        nlag = 1;
        mvar = nvar;
    case 'single'
        mvar = xvar+1;
end
ntrial = numel(x);
nlambda = numel(lambda);

% Compute covariance matrices
if arg.fast
    [CXX,CXY,XLAG] = olscovmat(x,y,lags,arg.type,arg.split,arg.zeropad,0);
else
    [CXX,CXY] = olscovmat(x,y,lags,arg.type,arg.split,arg.zeropad,1);
end

% Set up sparse regularization matrix
M = sparse(eye(mvar));
switch arg.method
    case 'ridge'
        M(1,1) = 0;
    case 'Tikhonov'
        M = M - 0.5*(diag(ones(mvar-1,1),1)+diag(ones(mvar-1,1),-1));
        M([mvar+2,end]) = 0.5;
        M([1,2,mvar+1]) = 0;
    case 'ols'
        lambda(:) = 0;
end
M = M/delta;

% Leave-one-out cross-validation
m = 0;
r = zeros(ntrial*arg.split,nlambda,yvar,nlag);
p = zeros(ntrial*arg.split,nlambda,yvar,nlag);
rmse = zeros(ntrial*arg.split,nlambda,yvar,nlag);
for i = 1:ntrial
    
    % Get segment size
    nseg = ceil(xobs(i)/arg.split);
    
    for j = 1:arg.split
        
        m = m+1;
        
        % Get segment indices
        iseg = nseg*(j-1)+1:min(nseg*j,xobs(i));
        
        if arg.fast % fast CV method
            
            % Validation set
            xlag = XLAG{m};
            
            % Training set
            Cxx = 0; Cxy = 0;
            itrain = 1:ntrial*arg.split;
            itrain(m) = [];
            for k = itrain
                Cxx = Cxx + CXX{k};
                Cxy = Cxy + CXY{k};
            end
            
        else % memory-efficient CV method
            
            % Validation set
            [xlag,ilag] = lagGen(x{i}(iseg,:),lags,arg.zeropad);
            xlag = [ones(numel(ilag),1),xlag]; %#ok<*AGROW>
            
            % Training set
            Cxx = CXX - xlag'*xlag;
            Cxy = CXY - xlag'*y{i}(iseg(ilag),:);
            
        end
        
        % Remove zero-padded indices
        if ~arg.zeropad
            iseg = iseg(1+max(0,lags(end)):end+min(0,lags(1)));
        end
        
        for n = 1:nlambda
            
            switch arg.type
                
                case 'multi'
                    
                    % Fit model
                    w = (Cxx + lambda(n)*M)\Cxy;
                    
                    % Predict output
                    pred = xlag*w;
                    
                    % Measure performance
                    [rt,pt] = corr(y{i}(iseg,:),pred);
                    r(m,n,:) = diag(rt);
                    p(m,n,:) = diag(pt);
                    rmse(m,n,:) = sqrt(mean((y{i}(iseg,:)-pred).^2,1))';
                    
                case 'single'
                    
                    ii = 0;
                    for jj = 2:xvar:nvar
                        
                        ii = ii+1;
                        
                        % Fit model
                        w = (Cxx(:,:,ii) + lambda(n)*M)\Cxy(:,:,ii);
                        
                        % Predict output
                        idx = [1,jj:jj+xvar-1];
                        pred = xlag(:,idx)*w;
                        
                        % Measure performance
                        [rt,pt] = corr(y{i}(iseg,:),pred);
                        r(m,n,:,ii) = diag(rt);
                        p(m,n,:,ii) = diag(pt);
                        rmse(m,n,:,ii) = sqrt(mean((y{i}(iseg,:) - ...
                            pred).^2,1))';
                        
                    end
                    
            end
            
        end
        
    end
    
end

% Format output arguments
if nargout > 3
    t = lags/fs*1e3;
end