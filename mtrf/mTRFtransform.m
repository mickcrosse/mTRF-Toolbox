function fmodel = mTRFtransform(bmodel,resp,varargin)
%MTRFTRANSFORM  Transform a backward model into a forward model.
%   FMODEL = MTRFTRANSFORM(BMODEL,RESP) transforms a backward decoding
%   model (stimulus to neural response) into the corresponding forward
%   encoding model (neural response to stimulus) as per Haufe et al.
%   (2014), enabling the neurophysiological interpretation of the backward
%   model weights. RESP are the neural responses that were used to compute
%   the original backward model.
%
%   MTRFTRANSFORM returns the transformed model in a structure with the
%   following fields:
%       'w'         -- transformed model weights (xvar-by-nlag-by-yvar)
%       't'         -- time lags (ms)
%       'fs'        -- sample rate (Hz)
%       'Dir'       -- direction of causality (forward=1, backward=-1)
%       'type'      -- type of model (multi-lag, single-lag)
%
%   If RESP is a matrix, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If RESP is a vector, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%
%   If RESP is a cell array containing multiple trials, the covariance
%   matrices of each trial are summed to produce one transformed model.
%
%   [...] = MTRFTRANSFORM(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default) or 2 to work
%                   along the rows.
%       'split'     A scalar specifying the number of segments in which to
%                   split each trial of data when computing the covariance
%                   matrices. This is useful for reducing memory usage on
%                   large datasets. By default, the entire trial is used.
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%       'verbose'   A numeric or logical specifying whether to execute in
%                   verbose mode: pass in 1 for verbose mode (default), or
%                   0 for non-verbose mode.
%
%   See also MTRFTRAIN, MTRFMULTITRAIN.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Haufe S, Meinecke F, Gorgen K, Dahne S, Haynes JD, Blankertz B,
%          Bieﬂmann F (2014) On the interpretation of weight vectors of
%          linear models in multivariate neuroimaging. NeuroImage
%          87:96-110.
%      [2] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <crossemj@tcd.ie>
%            Adam Bednar <bednara@tcd.ie>
%            Emily Teoh <teohe@tcd.ie>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Format data in cell arrays
[resp,nobs] = formatcells(resp,arg.dim,arg.split);
nfold = numel(resp);

% Convert time lags to samples
tmin = bmodel.t(1)*bmodel.fs/1e3;
tmax = bmodel.t(end)*bmodel.fs/1e3;

% Predict output
pred = mTRFpredict([],resp,bmodel,'dim',arg.dim,'zeropad',arg.zeropad,...
    'verbose',0);

% Convert to cell array
if ~iscell(pred)
    pred = {pred};
end

% Truncate output
if ~arg.zeropad
    resp = truncate(resp,tmin,tmax,nobs);
end

% Verbose mode
if arg.verbose
    v = verbosemode([],0,nfold);
end

% Initialize variables
Cxx = 0;
Css = 0;

for i = 1:nfold
    
    % Compute covariance matrices
    Cxx = Cxx + resp{i}'*resp{i};
    Css = Css + pred{i}'*pred{i};
    
    % Verbose mode
    if arg.verbose
        v = verbosemode(v,i,nfold);
    end
    
end

% Transform backward model weights
for i = 1:length(Css)
    bmodel.w(:,:,i) = Cxx*bmodel.w(:,:,i)/Css(i,i);
end
bmodel.w = permute(bmodel.w,[3,2,1]);

% Format output arguments
fmodel = struct('w',fliplr(bmodel.w),'t',-fliplr(bmodel.t),...
    'fs',bmodel.fs,'Dir',-bmodel.Dir,'type',bmodel.type);

% Verbose mode
if arg.verbose
    verbosemode(v,nfold+1,nfold,fmodel);
end

function v = verbosemode(v,fold,nfold,model)
%VERBOSEMODE  Execute verbose mode.
%   V = VERBOSEMODE(V,FOLD,NFOLD,MODEL) prints details about the progress
%   of the main function to the screen.

if fold == 0
    v = struct('msg',[],'h',[],'tocs',0);
    fprintf('\nTransform decoding model\n')
    fprintf('Computing covariance matrices\n')
    v.msg = ['%d/%d [',repmat(' ',1,nfold),']\n'];
    v.h = fprintf(v.msg,fold,nfold);
elseif fold <= nfold
    if fold == 1 && toc < 0.1
        pause(0.1)
    end
    v.tocs = v.tocs + toc;
    fprintf(repmat('\b',1,v.h))
    v.msg = ['%d/%d [',repmat('=',1,fold),repmat(' ',1,nfold-fold),'] - ',...
        '%.3fs/fold\n'];
    v.h = fprintf(v.msg,fold,nfold,v.tocs/fold);
    tic
end
if fold == nfold
    fprintf('Transforming model');
elseif fold > nfold
    fprintf(' - %.3fs\n',toc)
    modelsummary(model)
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

% Split data
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'split',1,validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'zeropad',true,validFcn); % zero-pad design matrix
addParameter(p,'verbose',true,validFcn); % verbose mode

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;