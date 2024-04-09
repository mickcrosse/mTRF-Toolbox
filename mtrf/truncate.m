function [y,yobs] = truncate(x,tmin,tmax,xobs)
%TRUNCATE  Truncate multivariate data in a cell array.
%   Y = TRUNCATE(X,TMIN,TMAX) returns the truncated version of X based on
%   the minimum and maximum time lags.
%
%   [Y,YOBS] = TRUNCATE(X,TMIN,TMAX) returns the number of observations
%   that were retained in each cell.
%
%   [...] = TRUNCATE(X,TMIN,TMAX,XOBS) specifies the number of observations
%   in each cell of X.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Set default values
if nargin < 4 || ~isempty(xobs)
    [~,xobs] = formatcells(x,1,1,0);
end

% Initialize variables
if ~iscell(x) && ~isempty(x)
    x = {x};
end
y = cell(size(x));
if nargout > 1
    yobs = zeros(size(x));
end

for i = 1:numel(x)
    
    % Get indices
    idx = max(0,tmax)+1:min(0,tmin)+xobs(i);
    
    % Truncate cell
    y{i} = x{i}(idx,:);
    
    if nargout > 1 % get dimensions
        yobs(i) = numel(idx);
    end
    
end