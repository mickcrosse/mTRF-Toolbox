function [y,rows,cols] = formatcells(x,dim,split,nanflag)
%FORMATCELLS  Format data in a cell array.
%   Y = FORMATCELLS(X,DIM) returns a cell array containing the data in X
%   with DIM as the first dimension. Pass in 2 for DIM to transpose the
%   data in each cell, or 1 to keep the existing dimensions. If X is a
%   vector or a matrix, it is returned within a 1-by-1 cell array. If X
%   contains row vectors, they are automatically transposed to column
%   vectors.
%
%   [Y,ROWS,COLS] = FORMATCELLS(...) returns the dimensions of the data in
%   each of the newly formatted cells.
%
%   [...] = FORMATCELLS(X,DIM,SPLIT) specifies the number of segments in
%   which to split the data in each cell. Values greater than 1 result in
%   in more cells containing smaller segments of data.
%
%   [...] = FORMATCELLS(X,DIM,SPLIT,NANFLAG) specifies whether to check the
%   data for NaN (Not-A-Number) values. Pass in 1 for NANFLAG to check for
%   NaNs (default), or 0 to skip check.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Set default values
if nargin < 2 || isempty(dim)
    dim = 1;
end
if nargin < 3 || isempty(split)
    split = 1;
end
if nargin < 4 || isempty(nanflag)
    nanflag = true;
end

% Initialize variables
if ~iscell(x) && ~isempty(x)
    x = {x};
end
ncells = numel(x)*split;
y = cell(ncells,1);
if nargout > 1
    rows = zeros(ncells,1);
end
if nargout > 2
    cols = zeros(ncells,1);
end
ii = 0;

for i = 1:numel(x)
    
    % Check for NaN values
    if nanflag
        if any(isnan(x{i}(:)))
            error('Input data contains NaN values.')
        end
    end
    
    % Orient data column-wise
    if dim == 2 || isrow(x{i})
        x{i} = x{i}';
    end
    
    % Compute max segment size
    if split > 1
        nobs = size(x{i},1);
        nseg = ceil(nobs/split);
    end
    
    for j = 1:split
        
        ii = ii+1;
        
        if split == 1 % use entire cell
            y{ii} = x{i};
        elseif split > 1 % split cell
            idx = nseg*(j-1)+1:min(nseg*j,nobs);
            y{ii} = x{i}(idx,:);
        end
        
        if nargout > 1 % get rows
            rows(ii) = size(y{ii},1);
        end
        if nargout > 2 % get columns
            cols(ii) = size(y{ii},2);
        end
        
    end
    
end