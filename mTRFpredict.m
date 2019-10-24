function [yhat,r,p,rmse] = mTRFpredict(stim,resp,w,fs,dir,tmin,tmax,b,varargin)
%MTRFPREDICT mTRF Toolbox model prediction.
%   YHAT = MTRFPREDICT(STIM,RESP,W,FS,DIR,TMIN,TMAX,B) returns the values
%   YHAT predicted by convolving the model weights W with the stimulus
%   features STIM or the brain responses RESP. To calculate the predictions
%   of an encoding model (stimulus to brain) pass in DIR=1, or to calculate
%   the predictions of a decoding model (brain to stimulus) pass in DIR=-1.
%   FS is a scalar containing the sample rate in Hertz. TMIN and TMAX are
%   scalars containing the minimum and maximum time lags in milliseconds
%   for generating time-lagged input features and should be the same as
%   those used to generate W. B is the bias term associated with W.
%
%   If STIM or RESP are multivariate, it is assumed that rows correspond to
%   observations and columns to variables. If they are univariate, it is
%   assumed that the first non-singleton dimension corresponds to
%   observations. STIM and RESP must have the same number of observations.
%
%   [...,R] = MTRFPREDICT(...) returns a scalar or a vector containing the
%   correlation coefficients between the predicted and observed values.
%
%   [...,P] = MTRFPREDICT(...) returns a scalar or a vector containing the
%   probabilities of the correlation coefficients.
%
%   [...,RMSE] = MTRFPREDICT(...) returns a scalar or a vector containing
%   the root-mean-square error between the predicted and observed values.
%
%   [...] = MTRFTRAIN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'dim'       a scalar specifying the dimension to work along: pass in 1
%               to work along the columns (default) or 2 to work along the
%               rows
%   'rows'      a string specifying the rows to use in the case of any
%               missing values (NaNs)
%                   'all'       use all rows, regardless of NaNs (default)
%                   'complete'  use only rows with no NaN values
%
%   See README for examples of use.
%
%   See also LAGGEN, MTRFTRAIN, MTRFTRANSFORM, MTRFCROSSVAL.
%
%   mTRF Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse, Giovanni Di Liberto, Edmund Lalor
%   Email: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 23-Oct-2019

% Decode input variable arguments
[dim,rows] = decode_varargin(varargin);

% Define X and Y
if tmin > tmax
    error('Value of TMIN must be < TMAX.')
elseif dir == 1
    x = stim; y = resp;
elseif dir == -1
    x = resp; y = stim;
    [tmin,tmax] = deal(tmax,tmin);
end
clear stim resp

% Arrange data column-wise
if dim == 2
    x = x'; y = y';
end
if size(x,1) == 1 && size(x,2) > 1
    x = x';
end
if size(y,1) == 1 && size(y,2) > 1
    y = y';
end
if size(x,1) ~= size(y,1)
    error('STIM and RESP must have the same number of observations.')
end

% Use only rows with no NaN values if specified
if strcmpi(rows,'complete')
    x = x(~any(isnan(y),2),:);
    y = y(~any(isnan(y),2),:);
    y = y(~any(isnan(x),2),:);
    x = x(~any(isnan(x),2),:);
elseif strcmpi(rows,'all') && (any(any(isnan(x))) || any(any(isnan(y))))
    error(['STIM or RESP missing values. Set argument ROWS to '...
        '''complete''.'])
end

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*dir);
tmax = ceil(tmax/1e3*fs*dir);

% Generate time-lagged features
X = [ones(size(x,1),1),lagGen(x,tmin:tmax)];

% Reformat and normalize model weights
w = [b;reshape(w,size(w,1)*size(w,2),size(w,3))]/fs;

% Compute prediction
yhat = X*w;

% Compute accuracy
if ~isempty(y)
    [r,p] = corr(y,yhat);
    r = diag(r); p = diag(p);
    rmse = sqrt(mean((y-yhat).^2,1))';
end

function [dim,rows] = decode_varargin(varargin)
%decode_varargin decode input variable arguments
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'dim')) && ~isempty(varargin{find(strcmpi(...
        varargin,'dim'))+1})
    dim = varargin{find(strcmpi(varargin,'dim'))+1};
    if ~isscalar(dim) || dim~=1 && dim~=2
        error(['Dimension argument must be a positive integer scalar '...
            'within indexing range.'])
    end
else
    dim = 1; % default: work along columns
end
if any(strcmpi(varargin,'rows')) && ~isempty(varargin{find(strcmpi(...
        varargin,'rows'))+1})
    rows = varargin{find(strcmpi(varargin,'rows'))+1};
    if ~any(strcmpi(rows,{'all','complete'}))
        error(['Invalid value for argument ROWS. Valid values are: '...
            '''all'', ''complete''.'])
    end
else
    rows = 'all'; % default: use all rows
end