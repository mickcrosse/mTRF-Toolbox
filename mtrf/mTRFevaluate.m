function [r,err] = mTRFevaluate(y,pred,varargin)
%MTRFEVALUATE  Evaluate the performance of a regression model.
%   R = MTRFEVALUATE(Y,PRED) returns the correlation between the predicted
%   output of a regression model PRED and the ground truth Y, based on
%   Pearson's linear correlation coefficient.
%
%   If Y or PRED are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   Y and PRED must have the same number of observations.
%
%   [R,ERR] = MTRFEVALUATE(Y,PRED) returns the error between the predicted
%   output and the ground truth, based on the mean squared error (MSE).
%
%   [...] = MTRFEVALUATE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both Y and PRED.
%       'corr'      A string specifying the correlation metric to use:
%                       'Pearson'   Pearson's linear correlation
%                                   coefficient (default): suitable for
%                                   data with a linear relationship
%                       'Spearman'  Spearman's rank correlation
%                                   coefficient: suitable for data with a
%                                   non-linear relationship
%       'error'     A string specifying the error metric to use:
%                       'mse'       mean squared error (default): take the
%                                   square root to convert to the original
%                                   units (i.e., RMSE)
%                       'mae'       mean absolute error: more robust to
%                                   outliers than MSE
%       'window'    A scalar specifying the window size over which to
%                   compute performance in samples. By default, the entire
%                   trial or segment is used.
%
%   See also CORR, CORRCOEF, TIEDRANK, IMMSE, MSE, MAE.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Orient data column-wise
if arg.dim == 2
    y = y';
    pred = pred';
end

% Get dimensions
[yobs,yvar] = size(y);
[pobs,pvar] = size(pred);
if pobs ~= yobs || pvar ~= yvar
    error(['Y and PRED arguments must have the same number of '...
        'observations and variables.'])
end
if arg.window
    nwin = floor(yobs/arg.window);
else
    nwin = 1;
end

% Initialize variables
r = zeros(nwin,yvar);
err = zeros(nwin,yvar);

for i = 1:nwin
    
    if arg.window % use window
        idx = arg.window*(i-1)+1:arg.window*i;
        yi = y(idx,:);
        predi = pred(idx,:);
        nobs = numel(idx);
    else % use entire trial
        yi = y;
        predi = pred;
        nobs = yobs;
    end
    
    % Compute error
    switch arg.error
        case 'mse'
            err(i,:) = sum(abs(yi-predi).^2,1)/nobs;
        case 'mae'
            err(i,:) = sum(abs(yi-predi),1)/nobs;
    end
    
    switch arg.corr
        case 'Spearman' % convert to rank values
            yi = num2rank(yi,nobs,yvar);
            predi = num2rank(predi,nobs,pvar);
    end
    
    % Demean signals
    y0 = bsxfun(@minus,yi,sum(yi,1)/nobs);
    pred0 = bsxfun(@minus,predi,sum(predi,1)/nobs);
    
    % Compute correlation coefficient
    r(i,:) = sum(y0.*pred0,1)./sqrt(sum(y0.^2,1).*sum(pred0.^2,1));
    
end

function xranked = num2rank(x,nobs,nvar)
%NUM2RANK  Rank numbers and average ties.
%   XRANKED = NUM2RANK(X) ranks the values in each column of X and averages
%   any tied ranks.

% Get dimensions
if nargin < 2
    nobs = size(x,1);
    nvar = size(x,2);
end
ranks = (1:nobs)';

% Initialize variables
xranked = zeros(nobs,nvar);

for i = 1:nvar
    
    % Sort data in ascending order
    [xsort,order] = sort(x(:,i));
    ranki = ranks;
    
    % Find ties
    ties = xsort(1:nobs-1) >= xsort(2:nobs);
    idx = [find(ties);nobs+2];
    maxt = numel(idx);
    
    % Average ties
    ctr = 1;
    while ctr < maxt
        m = idx(ctr); n = 2;
        while idx(ctr+1) == idx(ctr)+1
            ctr = ctr+1; n = n+1;
        end
        ranki(m:m+n-1) = sum(ranki(m:m+n-1))/n;
        ctr = ctr+1;
    end
    
    % Order ranks
    xranked(order,i) = ranki;
    
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

% Correlation metric
corrOptions = {'Pearson','Spearman'};
validFcn = @(x) any(validatestring(x,corrOptions));
addParameter(p,'corr','Pearson',validFcn);

% Error metric
errOptions = {'mse','mae'};
validFcn = @(x) any(validatestring(x,errOptions));
addParameter(p,'error','mse',validFcn);

% Window size
errorMsg = 'It must be a positive numeric scalar within indexing range.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'window',0,validFcn);

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially matched strings
arg.corr = validatestring(arg.corr,corrOptions);
arg.error = validatestring(arg.error,errOptions);