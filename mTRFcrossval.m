function [r,p,rmse,model,t] = mTRFcrossval(stim,resp,fs,map,tmin,tmax,lambda,tlims,nfolds)
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
%   [...,PRED,MODEL] = MTRFCROSSVAL(...) also returns the predictions PRED
%   and the linear mapping functions MODEL.
%
%   Inputs:
%   stim   - set of stimuli [cell{1,trials}(time by features)]
%   resp   - set of neural responses [cell{1,trials}(time by channels)]
%   fs     - sampling frequency (Hz)
%   map    - mapping direction (forward==1, backward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   lambda - ridge parameter values
%   tlims  - (optional) (NEW, NZ, 2019) specifies range or indexes of times 
%      that should be included in the model training and testing. If specific 
%      indexes are desired, then they should be specified in each cell of 
%      tlims, where the number of cells equals the number of trials.
%      (default: all indexes are used)
%      (see usetinds.m for more information on specifying tlims)
%   nfolds - (optional) specify the number of cross-validation folds 
%      (default: 10)
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
%   See also lagGen.m, mTRFTrain.m, mTRFpredict.m, mTRFcrossval.m.

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

% If tlims isn't specified, use all indexes
if nargin<9, tlims = []; end

% If nfolds isn't specified, do 10-fold cross-validation
if nargin<10, nfolds = 10; end

% Parse any addition function inputs
%if ~isempty(varargin),
%    for n = 2:2:length(varargin),
%        eval([varargin{n-1} '=varargin{n};']);
%    end
%end

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

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*map);
tmax = ceil(tmax/1e3*fs*map);

% Set up regularisation
dim1 = size(x{1},2)*length(tmin:tmax)+1;
dim2 = size(y{1},2);
model = zeros(numel(x),numel(lambda),dim1,dim2);
% If the inputs are a one-column array, do first-order Tikhonov regularization
% (Wong et al, 2018)
if size(x{1},2) == 1
    d = 2*eye(dim1);d([dim1+2,end]) = 1;
    u = [zeros(dim1,1),eye(dim1,dim1-1)];
    l = [zeros(1,dim1);eye(dim1-1,dim1)];
    M = d-u-l; M(:,1) = 0; M(1,:) = 0;
else % otherwise, do traditional ridge
    M = eye(dim1,dim1); M(1,1) = 0;
end

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
        tinds = tlims{i};
    else
        tinds = usetinds(tlims,fs,minlen);
    end
    x{i} = x{i}(tinds,:);
    y{i} = y{i}(tinds,:);
end
fprintf('Completed in %.3f s\n',toc(std_tm));

% Randomly resample time points
tottmidx = sum(cellfun(@(n) size(n,1),x));
rndtmsmp = randperm(tottmidx);
% Determine the starting index for each fold, with respect to the randomized samples
foldidx = (0:nfolds-1)*floor(tottmidx/nfolds);
    % the last fold has at most tottmidx/nfolds samples more than the others

disp('Starting model fitting...');
% Make the variables that will store the error and correlation values for each lambda and each fold
rmse = NaN(nfolds,numel(lambda),dim2); % root mean squared error
r = NaN(nfolds,numel(lambda),dim2); % Pearson's correlation
p = NaN(nfolds,numel(lambda),dim2); % significance of Pearson's correlation
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
    model = NaN(size(x{1},2),dim2,length(lambda));
    for j = 1:length(lambda)
        model(:,:,j) = (xtx+lambda(j)*M)\xty;
    end
    % Get the model lags (in ms)
    t = (tmin:tmax)/fs*1e3;
end

end