function [r,p,rmse,model] = mTRFcrossvalold(stim,resp,fs,map,tmin,tmax,lambda,tlims,varargin)
%mTRFcrossval mTRF Toolbox cross-validation function.
%   [R,P,RMSE] = MTRFCROSSVAL(STIM,RESP,FS,MAP,TMIN,TMAX,LAMBDA) performs
%   leave-one-out cross-validation on the set of stimuli STIM and the
%   neural responses RESP for the range of ridge parameter values LAMBDA.
%   As a measure of performance, it returns the correlation coefficients R
%   between the predicted and original signals, the corresponding p-values
%   P and the root-mean-square errors RMSE. Pass in MAP==1 to map in the
%   forward direction or MAP==-1 to map backwards. The sampling frequency
%   FS should be defined in Hertz and the time lags should be set in
%   milliseconds between TMIN and TMAX.
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
%   tlims  - (NEW, NZ, 2019) specifies range of time that should be included 
%      in the model training and testing. If specific indexes are desired,
%      then they should be specified in each cell of tlims, where the
%      number of cells equals the number of trials.
%      (default: all indexes are used)
%
%   Outputs:
%   r      - correlation coefficients
%   p      - p-values of the correlations
%   rmse   - root-mean-square errors
%   XX pred   - prediction [MAP==1: cell{1,trials}(lambdas by time by chans),
%            MAP==-1: cell{1,trials}(lambdas by time by feats)]
%          (NZ -- I would not recommend saving trial-by-trial predictions
%          here. At the end I compute the model based on all the training
%          data. Users should then create a prediction based on new data.)
%   model  - linear mapping function (MAP==1: trials by lambdas by feats by
%            lags by chans, MAP==-1: trials by lambdas by chans by lags by
%            feats)
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRAIN MTRFPREDICT MTRFMULTICROSSVAL.

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse, Giovanni Di Liberto, Edmund Lalor
%   Email: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 4-Feb-2019

% Initial parameters
nfolds = 10; % number of CV folds (NZ, 2019)

if nargin<9, tlims = []; end

if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

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
if size(x{1},2) == 1
    d = 2*eye(dim1);d([dim1+2,end]) = 1;
    u = [zeros(dim1,1),eye(dim1,dim1-1)];
    l = [zeros(1,dim1);eye(dim1-1,dim1)];
    M = d-u-l; M(:,1) = 0; M(1,:) = 0;
else
    M = eye(dim1,dim1); M(1,1) = 0;
end

%%% NZ edit -- creating trial-by-trial design matrices
disp('Creating the design matrices...');
% Create the design matrices
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
fprintf('** Completed design matrices in %.3f s\n',toc(std_tm));

%%% NZ edit -- Randomly resample time points, in order to determine the time indexes for
%%% each fold
tottmidx = sum(cellfun(@(n) size(n,1),x));
rndtmsmp = randperm(tottmidx);
% Determine the starting index for each fold
foldidx = (0:nfolds-1)*floor(tottmidx/nfolds);
    % the last fold has at most nfolds-1 samples more than the others

% Training
% disp('Training');
% X = cell(1,numel(x));
% for i = 1:numel(x)
%     % Generate lag matrix
%     X{i} = [ones(size(x{i},1),1),lagGen(x{i},tmin:tmax)];
%     % Calculate model for each lambda value
%     for j = 1:length(lambda)
%         model(i,j,:,:) = (X{i}'*X{i}+lambda(j)*M)\(X{i}'*y{i});
%     end
% end

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
    if n==nfolds,
        tstinds = rndtmsmp(foldidx(n)+1:end); % grab the rest of the samples if it's the last fold
    else
        tstinds = rndtmsmp(foldidx(n)+1:foldidx(n+1));
    end
    trninds = setxor(rndtmsmp,tstinds); % use the rest of the samples for training

    % Get the testing data
    xtst = cell_to_time_samples(x,tstinds);
    ytst = cell_to_time_samples(y,tstinds);
    [xtx,xty] = compute_linreg_matrices(x,y,trninds);
    
    % Calculate model for each lambda value
    fprintf('Cross-validating to each lambda (%d iterations)',length(lambda));
    for j = 1:numel(lambda)
        model_fold = (xtx+lambda(j)*M)\xty; % compute model
        pred_fold = xtst*model_fold; % compute prediction
        [rtmp,ptmp] = corr(ytst,squeeze(pred_fold));
        r(n,j,:) = diag(rtmp); p(n,j,:) = diag(ptmp);
        rmse(n,j,:) = sqrt(mean((ytst-pred_fold).^2,1));
        fprintf('.');
    end 
    fprintf('\n');
    clear xtx xty xtst ytst
    fprintf('Completed %.3f s\n',toc(mdlfitting_tm));
end

if nargout>3, % if the model is one of the outputs
    % Compute the final models for all training data using each lambda parameter
    disp('Computing models on all training data...');
    [xtx,xty] = compute_linreg_matrices(x,y);
    model = NaN(size(x{1},2),length(lambda),dim2);
    for j = 1:length(lambda)
        model(:,j,:) = (xtx+lambda(j)*M)\xty;
    end
end

% Testing
% pred = cell(1,numel(x));
% r = zeros(numel(x),numel(lambda),dim2);
% p = zeros(numel(x),numel(lambda),dim2);
% rmse = zeros(numel(x),numel(lambda),dim2);
% for i = 1:numel(x)
%     pred{i} = zeros(numel(lambda),size(y{i},1),dim2);
%     % Define training trials
%     trials = 1:numel(x);
%     trials(i) = [];
%     % Perform cross-validation for each lambda value
%     for j = 1:numel(lambda)
%         % Calculate prediction
%         pred{i}(j,:,:) = X{i}*squeeze(mean(model(trials,j,:,:),1));
%         % Calculate accuracy
%         [rtmp,ptmp] = corr(y{i},squeeze(pred{i}(j,:,:)));
%         r(i,j,:) = diag(rtmp); p(i,j,:) = diag(ptmp);
%         rmse(i,j,:) = sqrt(mean((y{i}-squeeze(pred{i}(j,:,:))).^2,1));
%     end
% end

end