function fmodel = mTRFtransform(bmodel,resp,varargin)
%MTRFTRANSFORM  mTRF model transformation.
%   FMODEL = MTRFTRANSFORM(BMODEL,RESP) transforms a backward decoding
%   model (stimulus to neural response) into the corresponding forward
%   encoding model (neural response to stimulus) as described in  Haufe et
%   al. (2014), enabling the neurophysiological interpretation of the
%   weights of the backward model. RESP are the neural responses that were
%   used to compute the original backward model.
%
%   MTRFTRANSFORM returns the transformed model in a structure with the
%   following fields:
%       'w'         -- transformed model weights (xvar-by-lags-by-yvar)
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
%
%   See mTRFdemos for examples of use.
%
%   See also MTRFTRAIN, MTRFPREDICT, MTRFCROSSVAL, MTRFMULTICROSSVAL.
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

%   Authors: Adam Bednar, Emily Teoh, Giovanni Di Liberto, Mick Crosse
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Apr 2014; Last revision: 10-Jan-2020

% Parse input arguments
arg = parsevarargin(varargin);

% Format data in cell arrays
resp = formatcells(resp,arg.dim);

% Predict output
pred = mTRFpredict([],resp,bmodel,'dim',arg.dim);

% Compute covariance matrices
covx = 0; covy = 0;
for i = 1:numel(resp)
    covx = covx + resp{i}'*resp{i};
    covy = covy + pred{i}'*pred{i};
end

% Transform backward model weights
bmodel.w = covx*bmodel.w/covy;

% Format output arguments
fmodel = struct('w',fliplr(bmodel.w),'t',-fliplr(bmodel.t),...
    'fs',bmodel.fs,'dir',-bmodel.dir,'type',bmodel.type);