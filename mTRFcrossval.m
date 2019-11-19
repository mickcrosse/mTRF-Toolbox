function [r,p,rmse,w,t,b] = mTRFcrossval(stim,resp,fs,map,tmin,tmax,lambda,varargin)
%mTRFcrossval mTRF Toolbox cross-validation function.
%   [R,P,RMSE] = MTRFCROSSVAL(STIM,RESP,FS,MAP,TMIN,TMAX,LAMBDA) performs
%   10-fold cross-validation on the set of stimuli STIM and the
%   neural responses RESP for the range of ridge parameter values LAMBDA.
%   As a measure of performance, it returns the correlation coefficients R
%   between the predicted and original signals, the corresponding p-values
%   P and the root-mean-square errors RMSE. Pass in MAP==1 to map in the
%   forward direction or MAP==-1 to map backwards. The sampling frequency
%   FS should be defined in Hertz and the time lags should be set in
%   milliseconds between TMIN and TMAX.
%     ** Important note **
%     Because each fold contains a random sampling of data, the
%     prediction accuracies may slightly vary each time this function is run.
%     Use these prediction accuracies to identify the optimal lambda value.
%     DO NOT USE THESE PREDICTION ACCURACIES TO REPORT MODEL PERFORMANCE!
%     Once the optimal model is identified, please use mTRFpredict.m to
%     compute model performance.
%
%   [...,W,T,B] = MTRFCROSSVAL(...) also returns the model W, the bias term
%   B, and the corresponding time array T for the model
%
%   Inputs:
%   stim   - set of stimuli [cell{1,trials}(time by features)]
%   resp   - set of neural responses [cell{1,trials}(time by channels)]
%   fs     - sampling frequency (Hz)
%   map    - mapping direction (forward==1, backward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   lambda - ridge parameter values
%   tlims (optional) - (NEW, NZ, 2019) specifies range or indexes of times 
%      that should be included in the model training and testing. If specific 
%      indexes are desired, then they should be specified in each cell of 
%      tlims, where the number of cells equals the number of trials.
%      Otherwise, set tlims=[] to use all of the data.
%      (default: all indexes are used)
%      (see usetinds.m for more information on specifying tlims)
%   nfolds - (optional) specify the number of cross-validation folds 
%      (default: 10)
%
%   Optional parameters (specify as 'parameter',value in function input)
%   Parameter   Value
%   'method'    a string specifying the regularization method to be used
%                   'Ridge'     Ridge regularization (default)
%                   'Tikhonov'  Tikhonov regularization
%   'scale'     a logical scalar specifying whether to scale regularization
%               parameter LAMBDA according to the data dimensions: pass in
%               0 to apply no scaling (default) or 1 to apply scaling
%   'dim'       a scalar specifying the dimension to work along: pass in 1
%               to work along the columns (default) or 2 to work along the
%               rows
%   'rows'      a string specifying the rows to use in the case of any
%               missing values (NaNs)
%                   'all'       use all rows, regardless of NaNs (default)
%                   'complete'  use only rows with no NaN values
%   'tlims'     specifies range or indexes of times that should be included 
%               in the model training and testing. If specific indexes are 
%               desired, then they should be specified in each cell of 
%               tlims, where the number of cells equals the number of trials.
%               Otherwise, set tlims=[] to use all of the data.
%               (default: all indexes are used)
%               (see usetinds.m for more information on specifying tlims)
%   'nfolds' -  specify the number of cross-validation folds 
%               (default: 10)
%
%   Outputs:
%   r      - correlation coefficients
%   p      - p-values of the correlations
%   rmse   - root-mean-square errors
%      ** Note: all prediction accuracy measures (r, p, rmse) are:
%                 MAP==1: fold by feats
%                 MAP==-1: fold by chans
%   model  - (optional) linear mapping function (MAP==1: lags by lambdas by feats, 
%            MAP==-1: lags by lambdas by chans)
%   t      - (optional) lags used by the model (in ms)
%
%   See README for examples of use.
%
%   See also lagGen.m, mTRFtrain.m, mTRFpredict.m, mTRFcrossval.m.

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2015) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse, Giovanni Di Liberto, Edmund Lalor
%   Recent edits (2019): Nathaniel J Zuk
%   Email: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 17-Oct-2019

% Decode input variable arguments
[method,dim,rows,nfolds,tlims] = decode_varargin(varargin);

% if the stim or response aren't cell arrays (only one trial was
% presented), make them a one element cell array
if ~iscell(stim), stim = {stim}; end
if ~iscell(resp), resp = {resp}; end

% Define x and y
if tmin > tmax
    error('Value of TMIN must be < TMAX')
end
if map == 1
    x = stim;
    y = resp;
elseif map == -1
    x = resp;
    y = stim;
    [tmin,tmax] = deal(tmax,tmin);
else
    error('Value of MAP must be 1 (forward) or -1 (backward)')
end
clear stim resp

size_chk = zeros(length(x),1); % flag if the # time samples in x and y aren't the same
for n = 1:length(x) % iterate through each trial
    % Arrange data column-wise
    if dim == 2
        x{n} = x{n}'; y{n} = y{n}';
    end
    % if it's a row array, flip to be a column array
    if size(x{n},1) == 1 && size(x{n},2) > 1
        x{n} = x{n}';
    end
    if size(y{n},1) == 1 && size(y{n},2) > 1
        y{n} = y{n}';
    end
    if size(x{n},1) ~= size(y{n},1)
%         error('Trial %d: STIM and RESP must have the same number of observations.',n)
        size_chk(n) = 1;
    end
end

% Warn the user if the # of time samples of any of the x-y pairs isn't the same
if sum(size_chk)>0
    warning(['STIM and RESP have a different number of time samples for at\n'...
        'least one of the trials. Arrays will be truncated to have the same\n'...
        'number of samples.']);
end

ncond = size(y{1},2); % number of conditions (columns of y) to fit

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*map);
tmax = ceil(tmax/1e3*fs*map);

ninputs = size(x{1},2); % get the number of columns in x, before it is 
    % replaced by the design matrix
disp('Creating the design matrices...');
% Create the design matrices trial by trial
std_tm = tic;
for i = 1:numel(x)
    % Generate lag matrix
    x{i} = [ones(size(x{i},1),1),lagGen(x{i},tmin:tmax)]; %%% skip constant term (NZ)
    % Set X and y to the same length
    minlen = min([size(x{i},1) size(y{i},1)]);
    x{i} = x{i}(1:minlen,:);
    y{i} = y{i}(1:minlen,:);
    % Remove time indexes, if specified
    if iscell(tlims), % if tlims is a cell array, it means that specific indexes were supplied
        tinds = usetinds(tlims{i},fs,minlen);
    else
        tinds = usetinds(tlims,fs,minlen);
    end
    x{i} = x{i}(tinds,:);
    y{i} = y{i}(tinds,:);
    % Use only rows with no NaN values if specified
    if strcmpi(rows,'complete')
        x{i} = x{i}(~any(isnan(y{i}),2),:);
        y{i} = y{i}(~any(isnan(y{i}),2),:);
        y{i} = y{i}(~any(isnan(x{i}),2),:);
        x{i} = x{i}(~any(isnan(x{i}),2),:);
    elseif strcmpi(rows,'all') && (any(any(isnan(x{i}))) || any(any(isnan(y{i}))))
        error(['STIM or RESP missing values. Set argument ROWS to '...
            '''complete''.'])
    end
end
fprintf('Completed in %.3f s\n',toc(std_tm));

% Create the regularization matrix
dim = size(x{1},2); % get the number of model parameters
if strcmpi(method,'Ridge') % Ridge regularization
    M = eye(dim,dim); M(1,1) = 0;
elseif strcmpi(method,'Tikhonov')  % Tikhonov regularization
    if size(x,2) > 1
        warning(['Tikhonov regularization may cause cross-channel '...
            'leakage for multivariate regression.'])
    end
    d = 2*eye(dim);d([dim+2,end]) = 1;
    u = [zeros(dim,1),eye(dim,dim-1)];
    l = [zeros(1,dim);eye(dim-1,dim)];
    M = d-u-l; M(:,1) = 0; M(1,:) = 0;
    M = M/2;
end

% Randomly resample time points
tottmidx = sum(cellfun(@(n) size(n,1),x));
rndtmsmp = randperm(tottmidx);
% Determine the starting index for each fold, with respect to the randomized samples
foldidx = (0:nfolds-1)*floor(tottmidx/nfolds);
    % the last fold has at most tottmidx/nfolds samples more than the others

disp('Starting model fitting...');
% Make the variables that will store the error and correlation values for each lambda and each fold
rmse = NaN(nfolds,numel(lambda),ncond); % root mean squared error
r = NaN(nfolds,numel(lambda),ncond); % Pearson's correlation
p = NaN(nfolds,numel(lambda),ncond); % significance of Pearson's correlation
% For each fold...
for n = 1:nfolds,
    mdlfitting_tm = tic;
    fprintf('Fold %d/%d\n',n,nfolds); % display which fold is being run
    % Identify the time samples used for this fold
    if n==nfolds, % if this is the last fold, get the rest of the indexes for testing
        tstinds = rndtmsmp(foldidx(n)+1:end); % grab the rest of the samples if it's the last fold
    else
        tstinds = rndtmsmp(foldidx(n)+1:foldidx(n+1));
    end
    trninds = setxor(rndtmsmp,tstinds); % use the rest of the samples for training

    % Get the testing data
    xtst = cell_to_time_samples(x,tstinds);
    ytst = cell_to_time_samples(y,tstinds);
    % Compute the matrices involved in the linear regression
    [xtx,xty] = compute_linreg_matrices(x,y,trninds);
    
    % Cross-validation
    fprintf('Cross-validating to each lambda (%d iterations)',length(lambda));
    for j = 1:numel(lambda)
        model_fold = (xtx+lambda(j)*M)\xty; % compute model
        pred_fold = xtst*model_fold; % compute prediction
        [rtmp,ptmp] = corr(ytst,squeeze(pred_fold)); % Pearson's
        r(n,j,:) = diag(rtmp); p(n,j,:) = diag(ptmp);
        rmse(n,j,:) = sqrt(mean((ytst-pred_fold).^2,1)); % root mean squared error
%         rmse(n,j,:) = std(ytst-pred_fold,1);
        fprintf('.');
    end 
    fprintf('\n');
    clear xtx xty xtst ytst
    fprintf('Completed in %.3f s\n',toc(mdlfitting_tm));
end

if nargout>3, % if the model is one of the outputs
    % Compute the final models for all training data using each lambda parameter
    disp('Computing model on all training data...');
    [xtx,xty] = compute_linreg_matrices(x,y);
    w = NaN(dim,ncond,length(lambda));
    for j = 1:length(lambda)
        w(:,:,j) = (xtx+lambda(j)*M)\xty;
    end
    w = w*fs;
    t = (tmin:tmax)/fs*1e3;
    b = w(1,:,:); % get the bias terms for each lambda
    w = reshape(w(2:end,:,:),[ninputs,length(t),ncond,length(lambda)]);
    % Get the model lags (in ms)
end

%% Decode varargin
function [method,dim,rows,nfolds,tlims] = decode_varargin(varargin)
%decode_varargin decode input variable arguments
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'method')) && ~isempty(varargin{find(strcmpi(...
        varargin,'method'))+1})
    method = varargin{find(strcmpi(varargin,'method'))+1};
    if ~any(strcmpi(method,{'Ridge','Tikhonov'}))
        error(['Invalid value for argument METHOD. Valid values are: '...
            '''Ridge'', ''Tikhonov''.'])
    end
else
    method = 'Ridge'; % default: use ridge method
end
if any(strcmpi(varargin,'dim')) && ~isempty(varargin{find(strcmpi(...
        varargin,'dim'))+1})
    dim = varargin{find(strcmpi(varargin,'dim'))+1};
    if ~isscalar(dim) || dim~=1 && dim~=2
        error(['Dimension argument must be a positive integer scalar '...
            'within indexing range.'])
    end
else
    dim = 1; % default: work along columns
end
if any(strcmpi(varargin,'rows')) && ~isempty(varargin{find(strcmpi(...
        varargin,'rows'))+1})
    rows = varargin{find(strcmpi(varargin,'rows'))+1};
    if ~any(strcmpi(rows,{'all','complete'}))
        error(['Invalid value for argument ROWS. Valid values are: '...
            '''all'', ''complete''.'])
    end
else
    rows = 'all'; % default: use all rows
end
% # folds
if any(strcmpi(varargin,'nfolds')) && ~isempty(varargin{find(strcmpi(...
        varargin,'nfolds'))+1})
    nfolds = varargin{find(strcmpi(varargin,'nfolds'))+1};
    if ~isinteger(nfolds)
        error('Invalid value for argument NFOLDS, it must be an integer')
    end
else
    nfolds = 10; % default: use all rows
end
% tlims
if any(strcmpi(varargin,'tlims')) && ~isempty(varargin{find(strcmpi(...
        varargin,'tlims'))+1})
    tlims = varargin{find(strcmpi(varargin,'tlims'))+1};
else
    tlims = []; % default: use all rows
end