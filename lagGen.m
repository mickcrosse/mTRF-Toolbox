function [xlag,idx] = lagGen(x,lags,zeropad,stack)
%LAGGEN  mTRF-Toolbox lag generator.
%   XLAG = LAGGEN(X,LAGS) returns a design matrix containing the time-
%   lagged features of input data X. X is a column vector or a matrix, with
%   the rows corresponding to observations and the columns to variables.
%   LAGS is a scalar or a vector of time lags in samples. Each lagged
%   feature is represented as a column in XLAG. If X is a matrix, each
%   column of XLAG is replaced with a matrix of variables, stacked in the
%   same order as they occur X.
%
%   [XLAG,IDX] = LAGGEN(X,LAGS) returns the indices of the design matrix
%   that were retained. If the design matrx is no zero-padded, this
%   information is useful for subsequent analysis.
%
%   XLAG = LAGGEN(X,LAGS,ZEROPAD) specifies whether to zero-pad the outer
%   rows of the design matrix or delete them. Pass in 1 to zero-pad them
%   (default), or 0 to delete them.
%
%   XLAG = LAGGEN(X,LAGS,ZEROPAD,STACK) specifies how to stack the lags and
%   variables in XLAG. Pass in 1 for STACK to stack variables in adjacent
%   columns for each lag (default), or 2 to stack lags in adjacent columns
%   for each variable. Note that the mTRF-Toolbox functions are designed to
%   unwrap model weights computed using option 1.
%
%   See also MTRFRESAMPLE, MTRFINTENSITY, MTRFPCA.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Apr 2014; Last revision: 10-Jan-2020

% Set default values
if nargin < 3 || isempty(zeropad)
    zeropad = true;
end
if nargin < 4 || isempty(stack)
    stack = 1;
end

% Get dimensions
nlag = numel(lags);
[nobs,nvar] = size(x);

% Generate time-lagged features
k = 1;
xlag = zeros(nobs,nvar*nlag);
for j = 1:nlag
    if stack == 1 % stack variables within lags
        if lags(j) < 0
            xlag(1:end+lags(j),k:k+nvar-1) = x(-lags(j)+1:end,:);
        elseif lags(j) > 0
            xlag(lags(j)+1:end,k:k+nvar-1) = x(1:end-lags(j),:);
        else
            xlag(:,k:k+nvar-1) = x;
        end
        k = k+nvar;
    elseif stack == 2 % stack lags within variables
        if lags(j) < 0
            xlag(1:end+lags(j),j:nlag:end) = x(-lags(j)+1:end,:);
        elseif lags(j) > 0
            xlag(lags(j)+1:end,j:nlag:end) = x(1:end-lags(j),:);
        else
            xlag(:,j:nlag:end) = x;
        end
    end
end

% Remove zero-padded rows
if zeropad
    idx = 1:nobs;
else
    idx = 1+max(0,lags(end)):nobs+min(0,lags(1));
    xlag = xlag(idx,:);
end