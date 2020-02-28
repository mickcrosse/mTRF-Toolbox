function [x,rows,cols] = formatcells(x,dim,nanflag)
%FORMATCELLS  Format data in cell array.
%   X = FORMATCELLS(X,DIM) returns a cell array containing the data in X
%   with DIM as the first dimension. Pass in 2 for DIM to transpose the
%   data in each cell, or 1 to keep the existing dimensions. If X is a
%   vector or a matrix, it is returned within a 1-by-1 cell array. If X
%   contains row vectors, they are automatically transposed to column
%   vectors.
%
%   [X,ROWS,COLS] = FORMATCELLS(...) returns the dimensions of the data in
%   each of the newly formatted cells.
%
%   [...] = FORMATCELLS(X,DIM,NANFLAG) specifies whether to check the data
%   for NaN (Not-A-Number) values. Pass in 1 for NANFLAG to check for NaNs
%   (default), or 0 to skip check.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Set default values
if nargin < 2
    dim = 1;
end
if nargin < 3
    nanflag = true;
end

if ~iscell(x) && ~isempty(x)
    x = {x};
end

if nargout > 1
    rows = zeros(size(x));
    cols = zeros(size(x));
end

for i = 1:numel(x)
    if nanflag
        if any(isnan(x{i}(:))) % check for NaN values
            error('Input data contains NaN values.')
        end
    end
    if dim == 2 || isrow(x{i}) % orient data column-wise
        x{i} = x{i}';
    end
    if nargout > 1 % get dimensions
        [rows(i),cols(i)] = size(x{i});
    end
end