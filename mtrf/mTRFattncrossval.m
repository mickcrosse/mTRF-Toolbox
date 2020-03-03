function [stats,stats1,stats2,t] = mTRFattncrossval(stim1,stim2,resp,fs,dir,tmin,tmax,lambda,varargin)
%MTRFATTNCROSSVAL  Cross-validation for attention decoder optimization.
%   STATS = MTRFATTNCROSSVAL(STIM1,STIM2,RESP,FS,DIR,TMIN,TMAX,LAMBDA)
%   cross validates a forward encoding model (stimulus to neural response)
%   or a backward decoding model (neural response to stimulus) over
%   multiple trials of data for optimizing an attention decoder. Models are
%   trained on the attended stimuli STIM1, and validated on both the
%   attended stimuli and unattended stimuli STIM2 as per O'Sullivan et al.
%   (2015). Pass in 1 for DIR to validate a forward model, or -1 to
%   validate a backward model. STIM1, STIM2 and RESP are cell arrays
%   containing corresponding trials of continuous data. FS is a scalar
%   specifying the sample rate in Hertz, and TMIN and TMAX are scalars
%   specifying the minimum and maximum time lags in milliseconds. For
%   backward models, MTRFATTNCROSSVAL automatically reverses the time lags.
%   LAMBDA is a scalar or vector of regularization values to be validated
%   and controls overfitting.
%
%   MTRFATTNCROSSVAL returns the cross-validation statistics for attention
%   decoding in a structure with the following fields:
%       'adi'       -- attention decoding index (ADI) based on the
%                      proportion of folds where the accuracy for attended
%                      stimuli was greater than that of the unattended
%                      stimuli (nlambda-by-yvar)
%       'dprime'    -- sensitivity index based on d', where the accuracy
%                      for attended stimuli is considered signal and that
%                      of the unattended stimuli is considered noise
%                      (nlambda-by-yvar)
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
%   validation statistics for the attended and unattended stimuli,
%   respectively, in structures with the following fields:
%       'acc'       -- prediction accuracy based on Pearson's correlation
%                      coefficient (ntrial-by-nlambda-by-yvar)
%       'err'       -- prediction error based on the mean squared error
%                      (ntrial-by-nlambda-by-yvar)
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
%       'type'      A string specifying type of model to fit:
%                       'multi'     use all lags simultaneously to fit a
%                                   multi-lag model (default)
%                       'single'    use each lag individually to fit
%                                   separate single-lag models
%       'acc'       A string specifying the accuracy metric to use:
%                       'Pearson'   Pearson's linear correlation
%                                   coefficient (default): suitable for
%                                   data with a linear relationship
%                       'Spearman'  Spearman's rank correlation
%                                   coefficient: suitable for data with a
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
%       'fast'      A numeric or logical specifying whether to use the fast
%                   cross-validation method (requires more memory) or the
%                   slower method (requires less memory): pass in 1 to use
%                   the fast method (default) or 0 to use the slower
%                   method. Note, both methods are numerically equivalent.
%
%   Notes:
%   Each iteration of MTRFATTNCROSSVAL partitions the N trials of data into
%   two subsets, fitting a model to N-1 trials (training set) and testing
%   on the left-out trial (validation set). Performance on the validation
%   set can be used to optimize hyperparameters (e.g., LAMBDA). Once
%   completed, it is recommended to test the model performance on separate
%   held-out data using the mTRFpredict function.
%
%   Discontinuous trials of data should not be concatenated prior to cross-
%   validation as this will introduce artifacts in places where time lags
%   cross over trial boundaries. Each trial should be input as an
%   individual cell of continuous data and MTRFATTNCROSSVAL will zero-pad
%   the trial boundaries appropriately.
%
%   See also CROSSVAL, DPRIME, MTRFPARTITION, MTRFTRAIN, MTRFPREDICT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] O'Sullivan JA, Power AJ, Mesgarani N, Rajaram S, Foxe JJ, Shinn-
%          Cunningham BG, Slaney M, Shamma SA, Lalor EC (2015) Attentional
%          Selection in a Cocktail Party Environment Can Be Decoded from
%          Single-Trial EEG. Cereb Cortex 25(7):1697-1706.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Edmund Lalor <edmundlalor@gmail.com>
%            Nate Zuk <zukn@tcd.ie>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Validate input parameters
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
    x = stim1; y = resp;
elseif dir == -1
    x = resp; y = stim1;
    [tmin,tmax] = deal(tmax,tmin);
else
    error('DIR argument must have a value of 1 or -1.')
end

% Format data in cell arrays
[x,xobs,xvar] = formatcells(x,arg.dim);
[y,yobs,yvar] = formatcells(y,arg.dim);
[z,zobs,~] = formatcells(stim2,arg.dim);

% Check equal number of observations
if ~isequal(xobs,yobs,zobs)
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
xvar = unique(xvar);
yvar = unique(yvar);
switch arg.type
    case 'multi'
        mvar = xvar*numel(lags)+1;
        nlag = 1;
    case 'single'
        mvar = xvar+1;
        nlag = numel(lags);
end
ntrial = numel(x);
nbatch = ntrial*arg.split;
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

% Initialize performance variables
acc1 = zeros(nbatch,nlambda,yvar,nlag);
err1 = zeros(nbatch,nlambda,yvar,nlag);
acc2 = zeros(nbatch,nlambda,yvar,nlag);
err2 = zeros(nbatch,nlambda,yvar,nlag);

% Leave-one-out cross-validation
n = 0;
for i = 1:ntrial
    
    % Max segment size
    seg = ceil(xobs(i)/arg.split);
    
    for j = 1:arg.split
        
        n = n+1;
        
        % Segment indices
        iseg = seg*(j-1)+1:min(seg*j,xobs(i));
        
        if arg.fast % fast method
            
            % Validation set
            xlag = XLAG{n};
            
            % Training set
            idx = 1:nbatch; idx(n) = [];
            Cxx = 0; Cxy = 0;
            for k = idx
                Cxx = Cxx + CXX{k};
                Cxy = Cxy + CXY{k};
            end
            
        else % memory-efficient method
            
            % Validation set
            [xlag,idx] = lagGen(x{i}(iseg,:),lags,arg.zeropad);
            xlag = [ones(numel(idx),1),xlag]; %#ok<*AGROW>
            
            % Training set
            Cxx = CXX - xlag'*xlag;
            Cxy = CXY - xlag'*y{i}(iseg(idx),:);
            
        end
        
        % Unattended validation set
        if dir == 1
            [zlag,idx] = lagGen(z{i}(iseg,:),lags,arg.zeropad);
            zlag = [ones(numel(idx),1),zlag];
        end
        
        % Remove zero-padded indices
        if ~arg.zeropad
            iseg = iseg(1+max(0,lags(end)):end+min(0,lags(1)));
        end
        
        for k = 1:nlambda
            
            switch arg.type
                
                case 'multi'
                    
                    % ---Attended Stimulus---
                    
                    % Fit linear model
                    w = (Cxx + lambda(k)*M)\Cxy;
                    
                    % Predict output
                    pred = xlag*w;
                    
                    % Evaluate performance
                    [acc1(n,k,:),err1(n,k,:)] = ...
                        mTRFevaluate(y{i}(iseg,:),pred,...
                        'acc',arg.acc,'err',arg.err);
                    
                    % ---Unattended Stimulus---
                    
                    if dir == 1
                        
                        % Predict output
                        pred = zlag*w;
                        
                        % Evaluate performance
                        [acc2(n,k,:),err2(n,k,:)] = ...
                            mTRFevaluate(y{i}(iseg,:),pred,...
                            'acc',arg.acc,'err',arg.err);
                        
                    elseif dir == -1
                        
                        % Evaluate performance
                        [acc2(n,k,:),err2(n,k,:)] = ...
                            mTRFevaluate(z{i}(iseg,:),pred,...
                            'acc',arg.acc,'err',arg.err);
                        
                    end
                    
                case 'single'
                    
                    for l = 1:nlag
                        
                        % ---Attended Stimulus---
                        
                        % Fit linear model
                        w = (Cxx(:,:,l) + lambda(k)*M)\Cxy(:,:,l);
                        
                        % Predict output
                        ilag = [1,xvar*(l-1)+2:xvar*l+1];
                        pred = xlag(:,ilag)*w;
                        
                        % Evaluate performance
                        [acc1(n,k,:,l),err1(n,k,:,l)] = ...
                            mTRFevaluate(y{i}(iseg,:),pred,...
                            'acc',arg.acc,'err',arg.err);
                        
                        % ---Unattended Stimulus---
                        
                        if dir == 1
                            
                            % Predict output
                            ilag = [1,xvar*(l-1)+2:xvar*l+1];
                            pred = zlag(:,ilag)*w;
                            
                            % Evaluate performance
                            [acc2(n,k,:,l),err2(n,k,:,l)] = ...
                                mTRFevaluate(y{i}(iseg,:),pred,...
                                'acc',arg.acc,'err',arg.err);
                            
                        elseif dir == -1
                            
                            % Evaluate performance
                            [acc2(n,k,:,l),err2(n,k,:,l)] = ...
                                mTRFevaluate(z{i}(iseg,:),pred,...
                                'acc',arg.acc,'err',arg.err);
                            
                        end
                        
                    end
                    
            end
            
        end
        
    end
    
end

% Compute ADI and d'
[adi,d] = mTRFattnevaluate(acc1,acc2);

% Format output arguments
stats = struct('adi',adi,'dprime',d);
stats1 = struct('acc',acc1,'err',err1);
stats2 = struct('acc',acc2,'err',err2);
if nargout > 2
    t = lags/fs*1e3;
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
addParameter(p,'fast',true,validFcn); % fast CV method
addParameter(p,'gpu',false,validFcn); % run on GPU

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially-matched strings
arg.method = validatestring(arg.method,regOptions);
arg.type = validatestring(arg.type,lagOptions);
arg.acc = validatestring(arg.acc,accOptions);
arg.err = validatestring(arg.err,errOptions);