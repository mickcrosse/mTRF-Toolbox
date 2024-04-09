function [xlag,idx] = lagGen(x,lags,zeropad,bias)
%LAGGEN  Generate time-lagged features of multivariate data.
%   XLAG = LAGGEN(X,LAGS) returns a design matrix containing the time-
%   lagged features of the input data X. X is a column vector or a matrix,
%   with the rows corresponding to observations and columns to variables.
%   LAGS is a scalar or a vector of time lags in samples. Each lagged
%   feature is represented as a column in XLAG. If X is a matrix, each
%   column of XLAG is replaced with a matrix of variables, stacked in the
%   same order as they occur X.
%
%   [XLAG,IDX] = LAGGEN(X,LAGS) returns the indices of the design matrix
%   that were retained. This information is useful if the design matrix is
%   not zero-padded (i.e., the end rows are deleted).
%
%   [...] = LAGGEN(X,LAGS,ZEROPAD) specifies whether to zero-pad the outer
%   rows of the design matrix or delete them. Pass in 1 for ZEROPAD to
%   zero-pad them (default), or 0 to delete them.
%
%   [...] = LAGGEN(X,LAGS,ZEROPAD,BIAS) specifies whether to include an
%   initial column of ones to allow for estimation of the bias term. Pass
%   in 1 for BIAS to include the bias term, or 0 to exclude it (default).
%
%   See also LAGMATRIX, MTRFRESAMPLE, MTRFENVELOPE.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Set default values
if nargin < 3 || isempty(zeropad)
    zeropad = true;
end
if nargin < 4 || isempty(bias)
    bias = false;
end

% Get dimensions
nlag = numel(lags);
[nobs,nvar] = size(x);

% Initialize variables
xlag = zeros(nobs,nvar*nlag);

% Generate time-lagged features
for i = 1:nlag
    ii = nvar*(i-1)+1:nvar*i;
    if lags(i) < 0
        xlag(1:end+lags(i),ii) = x(-lags(i)+1:end,:);
    elseif lags(i) > 0
        xlag(lags(i)+1:end,ii) = x(1:end-lags(i),:);
    else
        xlag(:,ii) = x;
    end
end

% Remove zero-padded rows
if zeropad
    idx = 1:nobs;
else
    idx = max(0,lags(end))+1:min(0,lags(1))+nobs;
    xlag = xlag(idx,:);
end

% Include bias term
if bias
    xlag = [ones(numel(idx),1),xlag];
end