function stats = mTRFnestedcrossval(stim,resp,fs,Dir,tmin,tmax,lambdas)
%MTRFNESTEDCROSSVAL Nested-loop leave-one-out cross-validation.
%   STATS = MTRFCROSSVAL(STIM,RESP,FS,DIR,TMIN,TMAX,LAMBDAS) runs a
%   nested-loop cross-validation. On each iteration, the model is trained
%   on all trials but one, and then tested on the left out trial.
%   MTRFCROSSVAL is used to identify the optimal lambda value in LAMBDAS
%   on the training set (not necessary if only one lambda value is
%   provided). Then, MTRFTRAIN is used to fit the model and MTRFPREDICT to
%   evaluate its predictive accuracy on the test set.
%
%   MTRFNESTEDCROSSVAL returns the cross-validation statistics in a 
%   structure with the following fields:
%       'acc'       -- prediction accuracy based on Pearson's correlation
%                      coefficient (nfold-by-nlambda-by-yvar)
%       'err'       -- prediction error based on the mean squared error
%                      (nfold-by-nlambda-by-yvar)
%
%   See also MTRFTRAIN, MTRFPREDICT, MTRFCROSSVAL.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Nate Zuk <zukn@tcd.ie>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% check if the inputs are cell arrays
if ~iscell(stim) || ~iscell(resp)
    error('STIM and RESP must be cell arrays containing more than one trial');
end

% check if there is more than one trial
if length(stim)<=1 || length(resp)<=1
    error('Both STIM and RESP must contain more than one cell');
end

% check if there are the same number of trials in stim and resp
if length(stim)~=length(resp)
    error('STIM and RESP must conatin the same number of cells');
end

% get the # of trials
nfold = numel(stim);

% get the # of conditions (different outputs vectors)
if Dir==1 
    ncond = size(resp{1},2);
elseif Dir==-1
    ncond = size(stim{1},2);
else
    error('MAP must be 1 (forward modeling) or -1 (backward modeling)');
end

% setup variables to store values on each iteration
r = NaN(nfold,ncond);
p = NaN(nfold,ncond);
rmse = NaN(nfold,ncond);
yhat = cell(nfold,1);
w = cell(nfold,1);
b = cell(nfold,1);
opt_lambda = NaN(nfold,1);

for n = 1:nfold
    
    fprintf('** Leaving out trial %d **\n',n)
    trtm = tic; % start the timer
    
    % Compute training and testing trials
    tst_tr = n;
    train_tr = setxor(1:nfold,n);

    % Do the model fitting
    if length(lambdas)==1 % if only one lambdas value is specified
        
        % use mTRFtrain
        model = mTRFtrain(stim(train_tr),resp(train_tr),fs,Dir,tmin,tmax,lambdas);
        
    else % otherwise, use mTRFcrossval
        
        [cvr,~,~,cvw,t,cvb] = mTRFcrossval(stim(train_tr),resp(train_tr),fs,Dir,tmin,tmax,lambdas,[]);
        stats = mTRFcrossval(stim(train_tr),resp(train_tr),fs,Dir,tmin,tmax,lambdas,[]);
        
        % compute the optimal lambdas
        tmp_r = squeeze(mean(mean(stats.r,3),1)); % average reconstruction accuracy across conditions and folds
        opt_idx = find(tmp_r==max(tmp_r),1,'first');
        
        % get the optimal model for testing
        opt_lambda(n) = lambdas(opt_idx);
        w{n} = cvw(:,:,:,opt_idx);
        b{n} = cvb(:,:,opt_idx);
        
    end
    
    % testing
    [~,stats(n)] = mTRFpredict(stim(tst_tr),resp(tst_tr),model);
    
    fprintf('** Completed training and testing in %.3f s\n',toc(trtm));
    
end