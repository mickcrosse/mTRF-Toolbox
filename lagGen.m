function xlag = lagGen(x,lags,stack)
%LAGGEN mTRF Toolbox lag generator.
%   XLAG = LAGGEN(X,LAGS) returns the matrix XLAG containing the time-
%   lagged features of X for the set of time lags in the vector LAGS. X is
%   a column vector or matrix, with rows corresponding to observations and
%   columns to variables. LAGS is a scalar or vector of positive and/or
%   negative integer values. Each lag feature is represented as a column of
%   XLAG. If X is multivariate, each column is replaced with a matrix of
%   variables, stacked in the same order as they occur X.
%
%   XLAG = LAGGEN(...,STACK) specifies how lags and variables are stacked
%   in XLAG. Pass in STACK=1 to stack variables in adjacent columns for
%   each lag (default) and pass in STACK=2 to stack lags in adjacent
%   columns for each variable. Note that MTRFTRAIN and MTRFPREDICT are
%   designed to unwrap model weights computed using option 1.
%
%   See also MTRFTRAIN, MTRFPREDICT, MTRFTRANSFORM, MTRFCROSSVAL.
%
%   mTRF Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 24-Oct-2019

% Set default feature order
if nargin < 3 || isempty(stack)
    stack = 1;
end

% Initialize variables
j = 1;
nfeat = size(x,2);
nlags = length(lags);
xlag = zeros(size(x,1),size(x,2)*length(lags));

% Generate time-lagged features
for i = 1:nlags
    if stack == 1 % stack variables within lags
        if lags(i) < 0
            xlag(1:end+lags(i),j:j+nfeat-1) = x(-lags(i)+1:end,:);
        elseif lags(i) > 0
            xlag(lags(i)+1:end,j:j+nfeat-1) = x(1:end-lags(i),:);
        else
            xlag(:,j:j+nfeat-1) = x;
        end
        j = j+nfeat;
    elseif stack == 2 % stack lags within variables
        if lags(i) < 0
            xlag(1:end+lags(i),i:nlags:end) = x(-lags(i)+1:end,:);
        elseif lags(i) > 0
            xlag(lags(i)+1:end,i:nlags:end) = x(1:end-lags(i),:);
        else
            xlag(:,i:nlags:end) = x;
        end
    else
        error(['Argument STACK must be an integer scalar of value 1 '...
            '(stack vars within lags) or 2 (stack lags within vars).'])
    end
end