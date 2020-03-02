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
%       'dir'       -- direction of causality (forward=1, backward=-1)
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
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%
%   See also MTRFTRAIN, MTRFMULTITRAIN.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.
%      [2] Haufe S, Meinecke F, Gorgen K, Dahne S, Haynes JD, Blankertz B,
%          Bieﬂmann F (2014) On the interpretation of weight vectors of
%          linear models in multivariate neuroimaging. NeuroImage
%          87:96-110.

%   Authors: Adam Bednar <bednara@tcd.ie>
%            Emily Teoh <teohe@tcd.ie>
%            Giovanni Di Liberto <diliberg@tcd.ie>
%            Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Format data in cell arrays
resp = formatcells(resp,arg.dim);

% Predict output
pred = mTRFpredict([],resp,bmodel,'dim',arg.dim,'zeropad',arg.zeropad);

% Convert to cell array
if ~iscell(pred)
    pred = {pred};
end

% Compute covariance matrices
Cxx = 0; Cyy = 0;
for i = 1:numel(resp)
    Cxx = Cxx + resp{i}'*resp{i};
    Cyy = Cyy + pred{i}'*pred{i};
end

% Transform backward model weights
bmodel.w = Cxx*bmodel.w/Cyy;

% Format output arguments
fmodel = struct('w',fliplr(bmodel.w),'t',-fliplr(bmodel.t),...
    'fs',bmodel.fs,'dir',-bmodel.dir,'type',bmodel.type);

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
addParameter(p,'zeropad',true,validFcn); % zero-pad design matrix
addParameter(p,'gpu',false,validFcn); % run on GPU

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;