function [r,p,mse,pred,model] = mTRFmulticrossval(stim,resp,resp1,resp2,fs,map,tmin,tmax,lambda1,lambda2)
%mTRFmulticrossval mTRF Toolbox multisensory cross-validation function.
%   [R,P,MSE] = MTRFMULTICROSSVAL(STIM,RESP,RESP1,RESP2,FS,MAP,TMIN,TMAX,
%   LAMBDA1,LAMBDA2) performs leave-one-out cross-validation of an
%   additive model for a multisensory dataset as follows:
%   1. Separate unisensory models are calculated using the set of stimuli
%      STIM and unisensory neural responses RESP1 and RESP2 for the range
%      of ridge parameter values LAMBDA1 and LAMBDA2 respectively.
%   2. The algebraic sums of the unisensory models for every combination of
%      LAMBDA1 and LAMBDA2 are calculated, i.e., the additive models.
%   3. The additive models are validated by testing them on the set of
%      multisensory neural responses RESP.
%   As a measure of performance, it returns the correlation coefficients R
%   between the predicted and original signals, the corresponding p-values
%   P and the mean squared errors MSE. The time lags T should be set in
%   milliseconds between TMIN and TMAX and the sampling frequency FS should
%   be defined in Hertz. Pass in MAP==1 to map in the forward direction or
%   MAP==-1 to map backwards. The neural responses in all three sensory
%   conditions must have been recorded for the same set of stimuli STIM.
%
%   [...,PRED,MODEL] = MTRFMULTICROSSVAL(...) also returns the predictions
%   PRED and the linear mapping functions MODEL.
%
%   Inputs:
%   stim    - set of stimuli [cell{1,trials}(time by features)]
%   resp    - set of multisensory neural responses [cell{1,trials}(time by channels)]
%   resp1   - set of unisensory 1 neural responses [cell{1,trials}(time by channels)]
%   resp2   - set of unisensory 2 neural responses [cell{1,trials}(time by channels)]
%   fs      - sampling frequency (Hz)
%   map     - mapping direction (forward==1, backward==-1)
%   tmin    - minimum time lag (ms)
%   tmax    - maximum time lag (ms)
%   lambda1 - unisensory 1 ridge parameter values
%   lambda2 - unisensory 2 ridge parameter values
%
%   Outputs:
%   r       - correlation coefficients
%   p       - p-values of the correlations
%   mse     - mean squared errors
%   pred    - prediction [MAP==1: cell{1,trials}(lambdas1 by lambdas2 by
%             time by chans), MAP==-1: cell{1,trials}(lambdas1 by lambdas2
%             by time by feats)]
%   model   - linear mapping function (MAP==1: trials by lambdas1 by
%             lambdas2 by feats by lags by chans, MAP==-1: trials by
%             lambdas1 by lambdas2 by chans by lags by feats)
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRAIN MTRFPREDICT MTRFCROSSVAL.

%   References:
%      [1] Crosse MC, Butler JS, Lalor EC (2015) Congruent visual speech
%          enhances cortical entrainment to continuous auditory speech in
%          noise-free conditions. J Neurosci 35(42):14195-14204.

%   Author: Michael Crosse
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Email: edmundlalor@gmail.com
%   Website: http://lalorlab.net/
%   April 2014; Last revision: 13 December 2016

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
dim1 = size(x{1},2)*length(tmin:tmax)+size(x{1},2);
dim2 = size(y{1},2);
model = zeros(numel(x),numel(lambda1),numel(lambda2),dim1,dim2);
if size(x{1},2) == 1
    d = 2*eye(dim1,dim1); d([1,end]) = 1;
    u = [zeros(dim1,1),eye(dim1,dim1-1)];
    l = [zeros(1,dim1);eye(dim1-1,dim1)];
    M = d-u-l;
else
    M = eye(dim1,dim1);
end

% Training
X = cell(1,numel(x));
for i = 1:numel(x)
    % Generate lag matrix
    X{i} = [ones(size(x{i})),lagGen(x{i},tmin:tmax)];
    if map == 1
        % Calculate unisensory models for each lambda value
        model1 = zeros(numel(lambda1),dim1,dim2);
        for j = 1:numel(lambda1)
            model1(j,:,:) = (X{i}'*X{i}+lambda1(j)*M)\X{i}'*resp1{i};
        end
        model2 = zeros(numel(lambda2),dim1,dim2);
        for j = 1:numel(lambda2)
            model2(j,:,:) = (X{i}'*X{i}+lambda2(j)*M)\X{i}'*resp2{i};
        end
    elseif map == -1
        % Generate lag matrices
        X1 = [ones(size(resp1{i})),lagGen(resp1{i},tmin:tmax)];
        X2 = [ones(size(resp2{i})),lagGen(resp2{i},tmin:tmax)];
        % Calculate unisensory models for each lambda value
        model1 = zeros(numel(lambda1),dim1,dim2);
        for j = 1:numel(lambda1)
            model1(j,:,:) = (X1'*X1+lambda1(j)*M)\X1'*y{i};
        end
        model2 = zeros(numel(lambda2),dim1,dim2);
        for j = 1:numel(lambda2)
            model2(j,:,:) = (X2'*X2+lambda2(j)*M)\X2'*y{i};
        end
        clear X1 X2
    end
    % Sum unisensory models for every combination of lambda values
    for j = 1:numel(lambda1)
        for k = 1:numel(lambda2)
            model(i,j,k,:,:) = model1(j,:,:)+model2(k,:,:);
        end
    end
    clear model1 model2
end
clear resp1 resp2

% Testing
pred = cell(1,numel(x));
r = zeros(numel(x),numel(lambda1),numel(lambda2),dim2);
p = zeros(numel(x),numel(lambda1),numel(lambda2),dim2);
mse = zeros(numel(x),numel(lambda1),numel(lambda2),dim2);
for i = 1:numel(x)
    pred{i} = zeros(numel(lambda1),numel(lambda2),size(y{i},1),dim2);
    % Define training trials
    trials = 1:numel(x);
    trials(i) = [];
    % Perform cross-validation for every combination of lambda values
    for j = 1:numel(lambda1)
        for k = 1:numel(lambda2)
            % Calculate prediction
            pred{i}(j,k,:,:) = X{i}*squeeze(mean(model(trials,j,k,:,:),1));
            % Calculate accuracy
            for l = 1:dim2
                [r(i,j,k,l),p(i,j,k,l)] = corr(y{i}(:,l),squeeze(pred{i}(j,k,:,l)));
                mse(i,j,k,l) = mean((y{i}(:,l)-squeeze(pred{i}(j,k,:,l))).^2);
            end
        end
    end
end

end