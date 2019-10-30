function [r,p,rmse,yhat,w,t,b,opt_lambda] = mTRFtrialbytrial(stim,resp,fs,map,tmin,tmax,lambda)
% Do trial-by-trial testing; on each iteration, the model is trained on all
% trials but one, and then tested on the left out trial. If lambda is just
% one value, then mTRFtrain is used to fit the model. Otherwise,
% mTRFcrossval is used, and the optimal model is selected that maximizes
% the cross-validation fit based on Pearson's correlation.

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
ntrials = length(stim);

% get the # of conditions (different outputs vectors)
if map==1 
    ncond = size(resp{1},2);
elseif map==-1
    ncond = size(stim{1},2);
else
    error('MAP must be 1 (forward modeling) or -1 (backward modeling)');
end

% setup variables to store values on each iteration
r = NaN(ntrials,ncond);
p = NaN(ntrials,ncond);
rmse = NaN(ntrials,ncond);
yhat = cell(ntrials,1);
w = cell(ntrials,1);
b = cell(ntrials,1);
opt_lambda = NaN(ntrials,1);
for n = 1:ntrials
    fprintf('** Leaving out trial %d **\n',n)
    trtm = tic; % start the timer
    
    % Compute training and testing trials
    tst_tr = n;
    train_tr = setxor(1:ntrials,n);

    % Do the model fitting
    if length(lambda)==1 % if only one lambda value is specified
        % use mTRFtrain
        [w{n},t,b{n}] = mTRFtrain(stim(train_tr),resp(train_tr),fs,map,tmin,tmax,lambda);
    else % otherwise, use mTRFcrossval
        [cvr,~,~,cvw,t,cvb] = mTRFcrossval(stim(train_tr),resp(train_tr),fs,map,tmin,tmax,lambda,[]);
        % compute the optimal lambda
        tmp_r = squeeze(mean(mean(cvr,3),1)); % average reconstruction accuracy across conditions and folds
        opt_idx = find(tmp_r==max(tmp_r),1,'first');
        % get the optimal model for testing
        opt_lambda(n) = lambda(opt_idx);
        w{n} = cvw(:,:,:,opt_idx);
        b{n} = cvb(:,:,opt_idx);
    end
    
    % testing
    [yhat{n},r(n,:),p(n,:),rmse(n,:)] = mTRFpredict(stim(tst_tr),resp(tst_tr),w{n},fs,map,tmin,tmax,b{n});
    
    fprintf('** Completed training and testing in %.3f s\n',toc(trtm));
end