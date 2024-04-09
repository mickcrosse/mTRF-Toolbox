function [y,t,cache] = mTRFresample(x,fsin,fsout,window,buff,varargin)
%MTRFRESAMPLE  Resample uniform data by a non-integer factor via averaging.
%   Y = MTRFRESAMPLE(X,FSIN,FSOUT) resamples the uniformly-sampled signal X
%   from a sample rate of FSIN to FSOUT by averaging the nearest neighbours
%   every FSIN/FSOUT samples. To avoid phase distortion, the inital window
%   is centered on the first sample. MTRFRESAMPLE does not apply an anti-
%   aliasing filter, but the data are in effect low-pass filtered by
%   convolution with a square window.
%
%   If X is a matrix, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If X is a vector, it is assumed that
%   the first non-singleton dimension corresponds to observations.
%
%   Y = MTRFRESAMPLE(X,FSIN,FSOUT,WINDOW) specifies the window size used to
%   average the nearest neighbours. Values greater than 1 result in overlap
%   between the data used to estimate neighbouring output frames, resulting
%   in increased smoothing. If WINDOW is not specified, it is set as 1 by
%   default. Note, changing the value of WINDOW does not affect the output
%   sample rate, only the degree of smoothing.
%
%   Y = MTRFRESAMPLE(X,FSIN,FSOUT,WINDOW,BUFF) concatenates a buffer of
%   initial data to the beginning of X to enable centering of the first
%   window at time t=0. The buffer should be passed from the final state of
%   previous data sampled at the input sample rate FSIN.
%
%   [Y,T] = MTRFRESAMPLE(...) returns the time axis of the newly resampled
%   signal Y.
%
%   [Y,T,CACHE] = MTRFRESAMPLE(...) caches the final state of the input
%   signal for initializing the buffer of subsequent data. This is useful
%   for resampling discrete segments of a continuous signal and real-time
%   applications.
%
%   [...] = MTRFRESAMPLE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default) or 2 to work
%                   along the rows.
%
%   See also RESAMPLE, DECIMATE, UPFIRDN, INTERP.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%            Edmund Lalor <edlalor@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Set default values
if nargin < 2 || isempty(fsin)
    fsin = 1;
end
if nargin < 3 || isempty(fsout)
    fsout = fsin;
end
if nargin < 4 || isempty(window)
    window = 1;
end
if nargin < 5
    buff = [];
end

% Arrange data column-wise
if arg.dim == 2 || isrow(x)
    x = x';
end

% Concatenate buffer
if ~isempty(buff)
    x = [buff;x];
    tau = size(buff,1);
else
    tau = 0;
end

% Get number of input/output frames
nin = size(x,1);
nout = round((nin-tau)/fsin*fsout);

% Create time axis
t = (0:nout-1)/fsout;

% Get half window size in seconds
n = 0.5*window/fsout;

% Compute new frames
y = zeros(nout,size(x,2));
if fsout ~= fsin || window > 1
    for i = 1:nout
        idx1 = max(0,round(fsin*(t(i)-n))+tau)+1;
        idx2 = min(nin,round(fsin*(t(i)+n))+tau+1);
        y(i,:) = mean(x(idx1:idx2,:),1);
    end
else
    x(1:tau) = [];
    y = x;
end

% Cache final state of input signal
cache = x(end-round(fsin*n)+1:end,:);

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

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;