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
%       'zeropad'   A numeric or logical specifying whether to zero-pad the
%                   outer rows of the design matrix or delete them: pass in
%                   1 to zero-pad them (default), or 0 to delete them.
%       'verbose'   A numeric or logical specifying whether to display
%                   details about transformation progress: pass in 1 to
%                   display details (default), or 0 to not display details.
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
nfold = numel(resp);

% Verbose mode
if arg.verbose
    fprintf('\nTransformation\n')
    fprintf('--------------\n')
end

% Predict output
pred = mTRFpredict([],resp,bmodel,'dim',arg.dim,'zeropad',arg.zeropad,...
    'verbose',0);

% Convert to cell array
if ~iscell(pred)
    pred = {pred};
end

% Verbose mode
if arg.verbose
    fprintf('Computing covariance matrices\n')
    msg = 'Fold %d/%d [';
    h = fprintf(msg,0,nfold);
    tocs = 0; tic
end

% Initialize variables
Cxx = 0;
Cyy = 0;

for i = 1:nfold
    
    % Compute covariance matrices
    Cxx = Cxx + resp{i}'*resp{i};
    Cyy = Cyy + pred{i}'*pred{i};
    
    % Verbose mode
    if arg.verbose
        fprintf(repmat('\b',1,h))
        msg = strcat(msg,'=');
        h = fprintf(msg,i,nfold);
        tocs = tocs+toc;
        if i == nfold
            fprintf('] - %.2fs/fold\n',tocs/nfold)
        else
            tic
        end
    end
    
end

% Verbose mode
if arg.verbose
    fprintf('Transforming model\n'); tic
end

% Transform backward model weights
bmodel.w = Cxx*bmodel.w/Cyy;

% Format output arguments
fmodel = struct('w',fliplr(bmodel.w),'t',-fliplr(bmodel.t),...
    'fs',bmodel.fs,'Dir',-bmodel.Dir,'type',bmodel.type);

% Verbose mode
if arg.verbose
    fprintf('model shape: %d x %d x %d - %.2fs\n',size(fmodel.w,1),...
        size(fmodel.w,2),size(fmodel.w,3),toc)
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
addParameter(p,'zeropad',true,validFcn); % zero-pad design matrix
addParameter(p,'verbose',true,validFcn); % print progress

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;