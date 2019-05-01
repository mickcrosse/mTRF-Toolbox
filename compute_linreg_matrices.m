function [xtx,xty] = compute_linreg_matrices(Xc,Yc,useidx,iteridxlim,varargin)
% Compute the autocovariance matrix xtx and covariance matrix xty for
% linear regression using the data Xc and Yc
% Inputs:
% - Xc = cell array of design matrixes X
% - Yc = cell array of output vectors (or matrixes) Y, corresponding to each cell in Xc
% - useidx = array of time indexes to use across all possible time indexes in Xc and Yc
%   (default: use all indexes)
% - iteridxlim = number of indexes to use on each iteration when computing
% the correlation matrixes (default: 10000)
% Outputs:
% - xtx, xty = corresponds to X'X and X'y used in linear regression
% Nate Zuk (2018)

verbose = 1;

% Put Xc and Yc in cells if they aren't cell array
if ~iscell(Xc), Xc = {Xc}; end
if ~iscell(Yc), Yc = {Yc}; end

% Use all indexes if not already specified
if nargin < 3 || isempty(useidx), useidx = 1:sum(cellfun(@(x) size(x,1),Xc)); end

% Use 10000 indexes on each iteration if not specified
if nargin < 4 || isempty(iteridxlim), iteridxlim = min([length(useidx) 10000]); end

if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

ndims = size(Xc{1},2); % number of dimensions in model
nouts = size(Yc{1},2); % number of output columns in each y

% Compute the X'X and X'y for linear regression
iteridx = [1:iteridxlim:length(useidx) length(useidx)+1]; % range of training indexes on each iteration
xtx = zeros(ndims);
xty = zeros(ndims,nouts);
if verbose, fprintf('Computing matrices in %d iterations',length(iteridx)-1); end
for ii = 2:length(iteridx),
    if verbose, fprintf('.'); end
    idx = useidx(iteridx(ii-1):iteridx(ii)-1); % training indexes
    x = cell_to_time_samples(Xc,idx); % compute design matrix for iteration
    y = cell_to_time_samples(Yc,idx); % compute output vector for iteration
    xtx = xtx + x'*x; % add to overall xtx
    xty = xty + x'*y; % add to overall xty
end  
if verbose, fprintf('\n'); end