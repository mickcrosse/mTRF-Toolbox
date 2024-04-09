function [y,t,cache] = mTRFenvelope(x,fsin,fsout,window,comp,buff,varargin)
%MTRFENVELOPE  Estimate the temporal envelope of an audio signal.
%   Y = MTRFENVELOPE(X,FSIN,FSOUT) computes the resampled temporal envelope
%   of the audio signal X by averaging the square of the nearest neighbours
%   every FSIN/FSOUT samples, taking the square root and logarithmically
%   scaling the RMS intensity (Lalor & Foxe, 2010).
%
%   If X is a matrix, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If it is a vector, it is assumed that
%   the first non-singleton dimension corresponds to observations.
%
%   Y = MTRFENVELOPE(X,FSIN,FSOUT,WINDOW) specifies the window size used to
%   average data. Values greater than 1 result in overlap between the data
%   used to estimate adjacent output frames resulting in increased envelope
%   smoothing. By default, a window size of 1 is used.
%
%   Y = MTRFENVELOPE(X,FSIN,FSOUT,WINDOW,COMP) specifies the amount of
%   compression applied to the envelope by raising the RMS value of X to 
%   the power of COMP. By default, a value of log10(2) is used to model
%   human auditory perception (Stevens, 1955).
%
%   Y = MTRFENVELOPE(X,FSIN,FSOUT,WINDOW,COMP,BUFF) concatenates a buffer
%   of initial data to the beginning of X to enable centering of the first
%   window at time t=0. The buffer should be passed from the final state of 
%   previous data sampled at the input sample rate FSIN.
%
%   [...] = MTRFENVELOPE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows.
%
%   See also ENVELOPE, HILBERT, RMS.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Lalor EC, Foxe JJ (2010) Neural responses to uninterrupted
%          natural speech can be extracted with precise temporal
%          resolution. Eur J Neurosci 31(1):189-193.
%      [2] Stevens SS (1955) The Measurement of Loudness. J Acoust Soc Am
%          27(2):815-829.

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
if nargin < 5 || isempty(comp)
    comp = log10(2);
end
if nargin < 6
    buff = [];
end

% Compute signal power
x = x.^2;

% Resample via moving average
[y,t,cache] = mTRFresample(x,fsin,fsout,window,buff,'dim',arg.dim);

% Apply compression to envelope
y = sqrt(y).^comp;

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