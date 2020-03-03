function [acc,err] = mTRFevaluate(y,pred,varargin)
%MTRFEVALUATE  Evaluate the prediction of a regression model.
%   ACC = MTRFEVALUATE(Y,PRED) evaluates the accuracy of the prediction
%   PRED relative to the ground truth Y based on Pearson's linear
%   correlation coefficient.
%
%   If Y or PRED are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   Y and PRED must have the same number of observations.
%
%   [ACC,ERR] = MTRFEVALUATE(Y,PRED) evaluates the error of the prediction
%   relative to the ground truth based on the mean squared error (MSE).
%
%   [...] = MTRFEVALUATE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both Y and PRED.
%       'acc'       A string specifying the accuracy metric to use:
%                       'Pearson'   Pearson's linear correlation
%                                   coefficient (default): suitable for
%                                   data with a linear relationship
%                       'Spearman'  Spearman's rank correlation
%                                   coefficient: suitable for data with a
%                                   non-linear relationship
%       'err'       A string specifying the error metric to use:
%                       'msc'       Mean squared error (default): take the
%                                   square root to convert to the original
%                                   units (i.e., RMSE)
%                       'mae'       Mean absolute error: more robust to
%                                   outliers than MSE
%
%   See also CORR, CORRCOEF, TIEDRANK, MSE, MAE.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Set default values
if nargin < 3 || isempty(arg.acc)
    arg.acc = 'Pearson';
end
if nargin < 4 || isempty(arg.err)
    arg.err = 'mse';
end

% Get dimensions
if arg.dim == 2
    y = y';
    pred = pred';
end
n = size(y,1);
if size(pred,1) ~= n
    error(['Y and PRED arguments must have the same number of '...
        'observations.'])
end

% Compute accuracy
switch arg.acc
    case 'Spearman'
        y = num2rank(y);
        pred = num2rank(pred);
end
muxy = sum(y).*sum(pred)/n;
sdxy = sqrt((sum(y.^2)-(sum(y).^2)/n).*(sum(pred.^2)-(sum(pred).^2)/n));
acc = (sum(y.*pred)-muxy)./sdxy;

% Compute error
switch arg.err
    case 'mse'
        exponent = 2;
    case 'mae'
        exponent = 1;
end
err = sum(abs(y-pred).^exponent,1)/n;

function xranked = num2rank(x)
%NUM2RANK  Rank numbers and average ties.
%   XRANKED = NUM2RANK(X) ranks the values in each column of X and averages
%   any tied ranks.

% Initialize variables
[nobs,nvar] = size(x);
xranked = zeros(nobs,nvar);

for i = 1:nvar
    
    % Sort data in ascending order
    [xsort,order] = sort(x(:,i));
    ranks = (1:nobs)';
    
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
        ranks(m:m+n-1) = sum(ranks(m:m+n-1))/n;
        ctr = ctr+1;
    end
    
    % Order ranks
    xranked(order) = ranks;
    
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

% Accuracy metric
accOptions = {'Pearson','Spearman'};
validFcn = @(x) any(validatestring(x,accOptions));
addParameter(p,'acc','Pearson',validFcn);

% Error metric
errOptions = {'mse','mae'};
validFcn = @(x) any(validatestring(x,errOptions));
addParameter(p,'err','msc',validFcn);

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially-matched strings
arg.acc = validatestring(arg.acc,accOptions);
arg.err = validatestring(arg.err,errOptions);