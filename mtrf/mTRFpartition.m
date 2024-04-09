function [strain,rtrain,stest,rtest] = mTRFpartition(stim,resp,k,testfold,varargin)
%MTRFPARTITION  Partition data into folds for cross-validation.
%   [STRAIN,RTRAIN] = MTRFPARTITION(STIM,RESP,K) partitions the stimulus
%   and response data into K equal folds for cross-validation. STIM and
%   RESP are vectors or matrices of continuous data and are returned as
%   K-by-1 cell arrays. To utilize all available data, the number of
%   samples in each fold is rounded up and adjusted for in the size of the
%   last fold. If K is not specified, it is set to 10 by default.
%
%   [STRAIN,RTRAIN,STEST,RTEST] = MTRFPARTITION(STIM,RESP,K,TESTFOLD) holds
%   out the fold specified by TESTFOLD and returns it as a separate test
%   set. If TESTFOLD is not specified, it is chosen at random by default.
%
%   [...] = MTRFPARTITION(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both STIM and RESP.
%       'equal'     A numeric or logical specifying whether to return folds
%                   of equal size or use all available data: pass in 1 for
%                   equal folds, or 0 to use all available data (default).
%
%   See also CVPARTITION, MTRFCROSSVAL.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Set default values
if nargin < 3 || isempty(k)
    k = 10;
end
if nargin < 4 || isempty(testfold)
    testfold = randi(k,1);
end

% Orient data column-wise
if arg.dim == 2
    stim = stim';
    resp = resp';
end

% Get dimensions
xobs = size(stim,1);
yobs = size(resp,1);

% Check equal number of observations
if ~isequal(xobs,yobs)
    error(['STIM and RESP arguments must have the same number of '...
        'observations.'])
end

% Define fold size
if arg.equal
    fold = floor(xobs/k);
else
    fold = ceil(xobs/k);
end

% Generate training set
strain = cell(k,1);
rtrain = cell(k,1);
for i = 1:k
    idx = fold*(i-1)+1:min(fold*i,xobs);
    strain{i} = stim(idx,:);
    rtrain{i} = resp(idx,:);
end

if nargout > 2
    
    % Generate test set
    stest = strain{testfold};
    rtest = rtrain{testfold};
    
    % Remove test set from training set
    strain(testfold) = [];
    rtrain(testfold) = [];
    
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

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'equal',false,validFcn); % equal fold size

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;