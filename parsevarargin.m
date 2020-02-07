function arg = parsevarargin(varargin)
%PARSEVARARGIN  Parse input arguments.
%   [PARAM1,PARAM2,...] = PARSEVARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   parses the input arguments of the main mTRF-Toolbox functions.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Jan 2020; Last revision: 05-Feb-2020

% Create parser object
p = inputParser;

% Dimension to work along
errorMsg = 'It must be a positive integer scalar within indexing range.';
validFcn = @(x) assert(x==1||x==2,errorMsg);
addParameter(p,'dim',1,validFcn);

% Regularization method
regOptions = {'ridge','Tikhonov','ols'};
validFcn = @(x) any(validatestring(x,regOptions));
addParameter(p,'method','ridge',validFcn);

% Model type
lagOptions = {'multi','single'};
validFcn = @(x) any(validatestring(x,lagOptions));
addParameter(p,'type','multi',validFcn);

% PCA algorithm
pcaOptions = {'eig','svd'};
validFcn = @(x) any(validatestring(x,pcaOptions));
addParameter(p,'pca','eig',validFcn);

% Segment data
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'split',1,validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'zeropad',true,validFcn); % zero pad design matrix
addParameter(p,'fast',true,validFcn); % fast cross-validation method
addParameter(p,'gpu',false,validFcn); % run on GPU
addParameter(p,'demean',true,validFcn); % demean data
addParameter(p,'plot',false,validFcn); % plot results
addParameter(p,'keepstate',false,validFcn); % keep state

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially-matched strings
arg.method = validatestring(arg.method,regOptions);
arg.type = validatestring(arg.type,lagOptions);
arg.pca = validatestring(arg.pca,pcaOptions);