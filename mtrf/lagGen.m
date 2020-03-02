function [xlag,idx] = lagGen(x,lags,zeropad,stack)
%LAGGEN  Generate time-lagged features of multivariate data.
%   XLAG = LAGGEN(X,LAGS) returns a design matrix containing the time-
%   lagged features of input data X. X is a column vector or a matrix, with
%   the rows corresponding to observations and the columns to variables.
%   LAGS is a scalar or a vector of time lags in samples. Each lagged
%   feature is represented as a column in XLAG. If X is a matrix, each
%   column of XLAG is replaced with a matrix of variables, stacked in the
%   same order as they occur X.
%
%   [XLAG,IDX] = LAGGEN(X,LAGS) returns the indices of the design matrix
%   that were retained. If the design matrix is not zero-padded (i.e., the
%   end rows are deleted), these data are useful for subsequent analysis.
%
%   XLAG = LAGGEN(X,LAGS,ZEROPAD) specifies whether to zero-pad the outer
%   rows of the design matrix or delete them. Pass in 1 to zero-pad them
%   (default), or 0 to delete them.
%
%   XLAG = LAGGEN(X,LAGS,ZEROPAD,STACK) specifies how to stack the lags and
%   variables in XLAG. Pass in 1 for STACK to stack variables in adjacent
%   columns for each lag (default), or 2 to stack lags in adjacent columns
%   for each variable. Note that mTRF-Toolbox functions are designed to
%   unwrap model weights computed using option 1.
%
%   See also LAGMATRIX, MTRFRESAMPLE, MTRFENVELOPE.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

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
xlag = zeros(nobs,nvar*nlag);
for i = 1:nlag
    if stack == 1 % stack variables within lags
        ilag = nvar*(i-1)+1:nvar*i;
        if lags(i) < 0
            xlag(1:end+lags(i),ilag) = x(-lags(i)+1:end,:);
        elseif lags(i) > 0
            xlag(lags(i)+1:end,ilag) = x(1:end-lags(i),:);
        else
            xlag(:,ilag) = x;
        end
    elseif stack == 2 % stack lags within variables
        if lags(i) < 0
            xlag(1:end+lags(i),i:nlag:end) = x(-lags(i)+1:end,:);
        elseif lags(i) > 0
            xlag(lags(i)+1:end,i:nlag:end) = x(1:end-lags(i),:);
        else
            xlag(:,i:nlag:end) = x;
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