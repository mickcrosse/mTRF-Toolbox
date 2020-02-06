function [pred,stats] = mTRFpredict(stim,resp,model,varargin)
%MTRFPREDICT  mTRF-Toolbox model prediction.
%   PRED = MTRFPREDICT(STIM,RESP,MODEL) predicts the output of a forward
%   encoding model (stimulus to neural response) or a backward decoding
%   model (neural response to stimulus) using time-lagged input features.
%   STIM and RESP are matrices or cell arrays containing corresponding
%   trials of continuous test data.
%
%   MODEL is a structure with the following fields:
%       'w'         -- normalized model weights (xvar-by-nlag-by-yvar)
%       'b'         -- normalized bias term (1-by-nlag-by-yvar)
%       't'         -- time lags (ms)
%       'fs'        -- sample rate (Hz)
%       'dir'       -- direction of causality (forward=1, backward=-1)
%       'type'      -- model type (multi-lag, single-lag)
%
%   If STIM or RESP are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   STIM and RESP must have the same number of observations.
%
%   If STIM or RESP are cell arrays containing multiple trials, separate
%   predictions are made for each trial and returned in a cell array. STIM
%   and RESP must contain the same number of trials.
%
%   [PRED,STATS] = MTRFPREDICT(...) returns a structure with the following
%   fields:
%       'r'         -- the correlation coefficients between the predicted
%                      and observed variables (ntrial-by-yvar-by-nlag)
%       'p'         -- the probabilities of the correlation coefficients
%       'rmse'      -- the root-mean-square error between the predicted and
%                      observed variables (ntrial-by-yvar-by-nlag)
%
%   [...] = MTRFPREDICT(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both STIM and RESP.
%       'split'     A scalar specifying the number of segments in which to
%                   split each trial of data when computing the covariance
%                   matrices. This is useful for reducing memory usage on
%                   large datasets. By default, the entire trial is used.
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%
%   See mTRFdemos for examples of use.
%
%   See also MTRFTRAIN, MTRFTRANSFORM, MTRFCROSSVAL, MTRFMULTICROSSVAL.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Lalor EC (2009) Modeling the human visual system using the
%          white-noise approach. 2009 4th International IEEE/EMBS
%          Conference on Neural Engineering, Antalya 589-592.

%   Authors: Mick Crosse, Giovanni Di Liberto, Edmund Lalor
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Apr 2014; Last revision: 05-Feb-2020

% Parse input arguments
arg = parsevarargin(varargin);

% Define X and Y variables
if model.dir == 1
    x = stim; y = resp;
elseif model.dir == -1
    x = resp; y = stim;
end

% Format data in cells column-wise
[x,xobs,xvar] = formatcells(x,arg.dim);
[y,yobs,yvar] = formatcells(y,arg.dim);
if yobs == 0 || yvar == 0
    yobs = xobs;
    yvar = size(model.w,3);
end

% Check equal number of observations
if ~isequal(xobs,yobs)
    error(['STIM and RESP arguments must have the same number of '...
        'observations.'])
end

% Convert time lags to samples
tmin = model.t(1)*model.fs/1e3;
tmax = model.t(end)*model.fs/1e3;
lags = tmin:tmax;

% Compute sampling interval
delta = 1/model.fs;

% Get dimensions
nlag = numel(lags);
xvar = unique(xvar);
yvar = unique(yvar);
nvar = xvar*nlag+1;
ntrial = numel(x);

% Format model weights
switch arg.type
    case 'multi'
        model.w = [model.b;reshape(model.w,[xvar*nlag,yvar])]*delta;
        nlag = 1;
    case 'single'
        model.w = [model.b;model.w]*delta;
end

% Test model
m = 0;
pred = cell(ntrial*arg.split,1);
r = zeros(ntrial*arg.split,yvar,nlag);
p = zeros(ntrial*arg.split,yvar,nlag);
rmse = zeros(ntrial*arg.split,yvar,nlag);
for i = 1:ntrial
    
    % Max segment size
    nseg = ceil(xobs(i)/arg.split);
    
    for j = 1:arg.split
        
        m = m+1;
        
        % Test data
        iseg = nseg*(j-1)+1:min(nseg*j,xobs(i));
        [xlag,ilag] = lagGen(x{i}(iseg,:),lags,arg.zeropad);
        xlag = [ones(numel(ilag),1),xlag]; %#ok<*AGROW>

        switch arg.type
            
            case 'multi'
                
                % Predict output
                pred{m} = xlag*model.w;
                
                % Measure performance
                if nargout > 1
                    [rt,pt] = corr(y{i}(ilag,:),pred{m});
                    r(m,:) = diag(rt);
                    p(m,:) = diag(pt);
                    rmse(m,:) = sqrt(mean((y{i}(ilag,:)-pred{m}).^2,1))';
                end
                
            case 'single'
                
                ii = 0;
                pred{m} = zeros(numel(ilag),yvar,nlag);
                for jj = 2:xvar:nvar
                    
                    ii = ii+1;
                    
                    % Predict output
                    idx = [1,jj:jj+xvar-1];
                    pred{m}(:,:,ii) = xlag(:,idx)*squeeze(model.w(:,ii,:));
                    
                    % Measure performance
                    if nargout > 1
                        [rt,pt] = corr(y{i}(ilag,:),pred{m}(:,:,ii));
                        r(m,:,ii) = diag(rt);
                        p(m,:,ii) = diag(pt);
                        rmse(m,:,ii) = sqrt(mean((y{i}(ilag,:)-pred{m}(:,:,ii)).^2,1))';
                    end
                    
                end
                
        end
                
    end
    
end

% Format output arguments
if ntrial == 1
    pred = pred{:};
end
if nargout > 1
    stats = struct('r',r,'p',p,'rmse',rmse);
end